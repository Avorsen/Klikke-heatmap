import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import tkinter as tk
from tkinter import filedialog, ttk, simpledialog
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
import scipy.interpolate
from matplotlib.path import Path
from scipy.spatial import ConvexHull
import os
import clipboard


click_handler = None
boundary_click_handler = None

class ClickHandler:
    def __init__(self, ax, text):
        self.ax = ax
        self.text = text
        self.cid = ax.figure.canvas.mpl_connect('button_press_event', self)

    def __call__(self, event):
        if event.inaxes != self.ax: return
        x = round(event.xdata, 2)  # round to two decimal places
        y = round(event.ydata, 2)
        z_value = simpledialog.askstring("Input", "Enter the Z value:")  # prompt user for z value

        if z_value is None:  # if user cancels the prompt
            return

        if z_value.startswith('-'):  # if z_value starts with '-'
            pass  # keep '-' in the beginning
        elif z_value.startswith('+'):  # if z_value starts with '+'
            pass  # keep '+' in the beginning
        else:  # if z_value does not start with '-' or '+'
            z_value = '-' + z_value  # add '-' to the beginning

        coordinates = [x, y, z_value]
        self.text.insert(tk.END, str(coordinates) + '\n')  # Insert new row with x, y, z
        self.ax.plot(x, y, 'ko')
        
        # Add the z_value as an annotation to the plot next to the point
        self.ax.annotate(z_value, (x, y), textcoords="offset points", xytext=(0,10), ha='center', color='red')
        
        self.ax.figure.canvas.draw()

class BoundaryClickHandler:
    def __init__(self, ax, coordinates):
        self.ax = ax
        self.coordinates = coordinates  # coordinates to choose boundary points from
        self.boundary = []
        self.cid = ax.figure.canvas.mpl_connect('button_press_event', self)
        self.boundary_closed = False

    def __call__(self, event):
        if self.boundary_closed:
            return

        if event.inaxes != self.ax:
            return

        x = round(event.xdata, 2)
        y = round(event.ydata, 2)

        closest_coord = min(self.coordinates, key=lambda coord: (coord[0] - x)**2 + (coord[1] - y)**2)

        if self.boundary and (self.boundary[0] == closest_coord):
            self.boundary_closed = True
            print("Boundary has been set.")
            # Draw line from last point to first point to close the loop
            self.ax.plot(*zip(self.boundary[-1], self.boundary[0]), color='r')
            # Add text to the top left corner
            self.ax.text(0.05, 0.95, 'Perimeter for plantegning satt', transform=self.ax.transAxes, color='red', fontsize=14, va='top')
            self.ax.text(0.05, 0.9, 'Du kan nå trykke på generer', transform=self.ax.transAxes, color='red', fontsize=14, va='top')
            self.ax.figure.canvas.draw()
            return

        self.boundary.append(closest_coord)
        self.ax.plot(*closest_coord, 'ro')  # plot boundary point

        if len(self.boundary) > 1:
            self.ax.plot(*zip(*self.boundary[-2:]), color='r')  # plot boundary line

        self.ax.figure.canvas.draw()

    
def set_boundary():
    coordinates = coordinates_text.get("1.0", tk.END).strip().split('\n')
    coordinates = [eval(coord) for coord in coordinates]
    
    xy_coordinates = [(x, y) for x, y, _ in coordinates]

    fig, ax = plt.subplots(figsize=(10,10))  # Define size of the figure here
    ax.scatter(*zip(*xy_coordinates), color='black')

    global boundary_click_handler
    boundary_click_handler = BoundaryClickHandler(ax, xy_coordinates)

    # Create a new top level window
    top = tk.Toplevel(root)
    top.title("Set Boundary")  # You can set title for the new window

    # Create the canvas on this new window
    canvas = FigureCanvasTkAgg(fig, master=top)
    canvas.draw()
    canvas.get_tk_widget().pack()  # Using pack() for simplicity




    def undo(self):
        if self.text.get('1.0', tk.END).strip():
            lines = self.text.get('1.0', tk.END).strip().split('\n')
            lines = lines[:-1]
            self.text.delete('1.0', tk.END)
            self.text.insert('1.0', '\n'.join(lines) + '\n')
            self.ax.lines[-1].remove()
            self.ax.figure.canvas.draw()

def upload_image():
    filepath = filedialog.askopenfilename()
    image = np.array(plt.imread(filepath))

    fig, ax = plt.subplots(figsize=(10,10))  # Define size of the figure here
    ax.imshow(image)
    ax.invert_yaxis()  # flip the y-axis

    global click_handler
    click_handler = ClickHandler(ax, coordinates_text)

    # Create a new top level window
    top = tk.Toplevel(root)
    top.title("Image Window")  # You can set title for the new window

    # Create the canvas on this new window
    canvas = FigureCanvasTkAgg(fig, master=top)
    canvas.draw()
    canvas.get_tk_widget().pack()  # Using pack() for simplicity

def reproduce_drawing():
    # Create a new figure for the reproduced points
    fig, ax = plt.subplots(figsize=(10,10))

    # Create a temporary ClickHandler for the reproduced points
    tmp_handler = ClickHandler(ax, coordinates_text)

    coordinates = coordinates_text.get("1.0", tk.END).strip().split('\n')
    for coordinate in coordinates:
        x, y, z = eval(coordinate)
        tmp_handler.ax.plot(x, y, 'ko')

    # Create a new top level window
    top = tk.Toplevel(root)
    top.title("Vis XY")  # You can set title for the new window

    # Create the canvas on this new window
    canvas = FigureCanvasTkAgg(fig, master=top)
    canvas.draw()
    canvas.get_tk_widget().pack()  # Using pack() for simplicity

    tmp_handler.ax.figure.canvas.draw()


# Create and register custom colormaps
colors = ["DarkOrange", "PaleGoldenrod", "white"]
cmap = mcolors.LinearSegmentedColormap.from_list("rosegold", colors)
plt.register_cmap(cmap=cmap)

colors2 = ["chocolate", "DarkOrange", "PaleGoldenrod", "white"]
cmap2 = mcolors.LinearSegmentedColormap.from_list("rosegold chocolate", colors2)
plt.register_cmap(cmap=cmap2)

colors3 = ["brown", "orange", "white"]
cmap3 = mcolors.LinearSegmentedColormap.from_list("brown_orange_white", colors3)
plt.register_cmap(cmap=cmap3)

colors4 = ["sienna", "orange", "white"]
cmap4 = mcolors.LinearSegmentedColormap.from_list("brownred_orange_white", colors4)
plt.register_cmap(cmap=cmap4)

colors5 = ["sienna", "yellow", "white"]
cmap5 = mcolors.LinearSegmentedColormap.from_list("brownred_yellow_white", colors5)
plt.register_cmap(cmap=cmap5)

colors6 = ["darkorange", "goldenrod", "palegoldenrod", "white"]
cmap6 = mcolors.LinearSegmentedColormap.from_list("goldenrod", colors6)
plt.register_cmap(cmap=cmap6)




# Create an opace of the Ylord
original_cmap = plt.get_cmap("YlOrRd_r")
colors = original_cmap(np.arange(original_cmap.N))
alpha = 0.8 #rediger denne for å sette opace styrken
colors[:, -1] = alpha  # Change alpha values of all colors
new_cmap_name = "YlOrRd_r_opace"
new_cmap = mcolors.LinearSegmentedColormap.from_list(new_cmap_name, colors)
plt.register_cmap(cmap=new_cmap)


# Create an opace version of the "rosegold chocolate"
original_cmap2 = plt.get_cmap("rosegold chocolate")
colors2_opace = original_cmap2(np.arange(original_cmap2.N))
alpha = 0.84 # Modify this to set opacity level
colors2_opace[:, -1] = alpha  # Change alpha values of all colors
new_cmap2_name = "rosegold chocolate_opace"
new_cmap2 = mcolors.LinearSegmentedColormap.from_list(new_cmap2_name, colors2_opace)
plt.register_cmap(cmap=new_cmap2)

def generate_heatmap():

    coordinates = coordinates_text.get("1.0", tk.END).strip().split('\n')
    points = []
    for coordinate in coordinates:
        x, y, z = eval(coordinate)
        points.append((x, y, int(z)))

    x, y, z = zip(*points)

    # Create grid
    xi = np.linspace(min(x), max(x), 1000)
    yi = np.linspace(min(y), max(y), 1000)
    xi, yi = np.meshgrid(xi, yi)

    # Interpolate
    zi = scipy.interpolate.griddata((x, y), z, (xi, yi), method='cubic')

    boundary_path = None
    points2D = np.array(points)[:,:2] # we only need x and y for 2D hull

    if boundary_click_handler and boundary_click_handler.boundary_closed:
        boundary_path = Path(boundary_click_handler.boundary)
    else:
        # Create convex hull
        hull = ConvexHull(points2D)
        boundary_path = Path(hull.points[hull.vertices])

    # Interpolate
    zi = scipy.interpolate.griddata((x, y), z, (xi, yi), method='cubic')

    # Create a mask for the region outside the contours
    points = np.column_stack((xi.flatten(), yi.flatten()))
    mask = ~boundary_path.contains_points(points).reshape(xi.shape)

    masked_zi = np.ma.array(zi, mask=mask)

    fig, ax = plt.subplots(figsize=(10,10))

    # Make the background transparent
    fig.patch.set_visible(False)
    ax.axis('off')

    min_colormap_range = float(min_colormap_range_entry.get()) if min_colormap_range_entry.get() else min(z) - 10
    max_colormap_range = float(max_colormap_range_entry.get()) if max_colormap_range_entry.get() else max(z)
    
    im = ax.imshow(masked_zi, extent=(min(x), max(x), min(y), max(y)), cmap=selected_colormap.get(), origin='lower', vmin=min_colormap_range, vmax=max_colormap_range)
    # Colorbar
    cbar = plt.colorbar(im, ax=ax, aspect=len(yi)//25, shrink=0.8)
    cbar.ax.set_ylabel('Setningsmålinger (mm)', rotation=-90, va="bottom")


    # Contour lines
    is_contour_lines_enabled = contour_lines_enabled.get()
    contour_line_width = float(contour_lines_width_entry.get())
    num_contour_lines = int(contour_lines_num_entry.get())

    if is_contour_lines_enabled:
        contour_interval = (max(z) - min(z)) / num_contour_lines
        contour_levels = np.arange(min(z), max(z), contour_interval)
        ax.contour(xi, yi, masked_zi, levels=contour_levels, colors='gray', linewidths=contour_line_width, linestyles='solid')



    # Save the figure without z values
    i = 1
    while True:
        saved_image_name = f"varmekart_{i}" if i > 1 else "varmekart"
        if not os.path.exists(os.path.join(os.path.expanduser("~"), "Desktop", f"{saved_image_name}.png")):
            break
        i += 1

    plt.savefig(os.path.join(os.path.expanduser("~"), "Desktop", f"{saved_image_name}.png"), dpi=300, transparent=True, bbox_inches='tight', pad_inches=0)

    # Display z values on the image
    font_size = 8  # or another size if you prefer
    for x_coord, y_coord, z_value in zip(x, y, z):
        z_value_str = str(int(z_value)) if float(z_value).is_integer() else str(z_value)
        ax.text(x_coord, y_coord, z_value_str, fontsize=font_size, ha='center', va='center', color='black')

    # Save the figure with z values
    plt.savefig(os.path.join(os.path.expanduser("~"), "Desktop", f"{saved_image_name}_with_z_values.png"), dpi=300, bbox_inches='tight', pad_inches=0)
    plt.show()

    print("Heatmap saved to Desktop. You can find the images with and without z-values.")




root = tk.Tk()


ascii_art = """
 _   _            _                          ______ _                  _ _  __ _   
| | | |          | |                         |  ___| |                | (_)/ _| |  
| |_| | ___  __ _| |_ _ __ ___   __ _ _ __   | |_  | | ___   ___  _ __| |_| |_| |_ 
|  _  |/ _ \/ _` | __| '_ ` _ \ / _` | '_ \  |  _| | |/ _ \ / _ \| '__| | |  _| __|
| | | |  __/ (_| | |_| | | | | | (_| | |_) | | |   | | (_) | (_) | |  | | | | | |_ 
\_| |_/\___|\__,_|\__|_| |_| |_|\__,_| .__/  \_|   |_|\___/ \___/|_|  |_|_|_|  \__|
                                     | |                                           
                                     |_|                                           
"""

ascii_art_label = tk.Label(root, text=ascii_art, font=("Courier", 12))
ascii_art_label.grid(row=0, column=0, columnspan=2, pady=(10, 10))


upload_button = tk.Button(root, text="Velg Plantegning", command=upload_image)
upload_button.grid(row=1, column=0)

boundary_button = tk.Button(root, text="Set Boundary", command=set_boundary)
boundary_button.grid(row=2, column=0) 

undo_button = tk.Button(root, text="Undo", command=lambda: click_handler.undo())
undo_button.grid(row=1, column=1)

coordinates_text = tk.Text(root, width=40, height=15)  # create a Text widget
coordinates_text.grid(row=3, column=0, columnspan=3, pady=20)

# Prepopulate the text box with house coordinates
house_coordinates = [
[129.25, 519.22, '0'],
[349.84, 508.4, '-35'],
[533.66, 497.59, '-53'],
[695.85, 495.43, '-70'],
[754.24, 402.44, '-75'],
[814.79, 279.17, '-81'],
[892.65, 181.85, '-89'],
[953.2, 110.49, '-95'],
[1149.99, 101.84, '-115'],
[1266.77, 84.54, '-127'],
[1390.04, 80.21, '-139'],
[1521.96, 84.54, '-152'],
[1571.7, 164.55, '-157'],
[1573.86, 283.5, '-157'],
[1569.53, 376.49, '-157'],
[1565.21, 484.62, '-157'],
[1573.86, 631.67, '-157'],
[1584.67, 824.14, '-158'],
[1576.02, 984.17, '-158'],
[1567.37, 1113.93, '-157'],
[1534.93, 1315.05, '-153'],
[1506.82, 1397.23, '-151'],
[1359.76, 1390.74, '-136'],
[1217.03, 1388.57, '-122'],
[1102.42, 1390.74, '-110'],
[970.5, 1397.23, '-97'],
[810.47, 1397.23, '-81'],
[752.08, 1388.57, '-75'],
[667.74, 1338.84, '-67'],
[693.69, 1260.98, '-69'],
[708.83, 1148.53, '-71'],
[708.83, 1010.12, '-71'],
[708.83, 936.6, '-71'],
[570.42, 956.06, '-57'],
[464.45, 951.73, '-46'],
[367.14, 945.25, '-37'],
[256.85, 921.46, '-26'],
[133.58, 960.38, '0'],
[133.58, 767.91, '0'],
[133.58, 594.91, '0'],
[354.16, 722.5, '-35'],
[544.47, 716.01, '-54'],
[691.53, 705.2, '-69'],
[864.53, 545.17, '-86'],
[907.78, 703.04, '-91'],
[877.51, 882.53, '-88'],
[875.35, 1020.94, '-88'],
[866.69, 1176.64, '-87'],
[875.35, 1276.12, '-88'],
[1104.58, 1273.96, '-110'],
[1271.1, 1267.47, '-127'],
[1383.55, 1265.31, '-138'],
[1372.74, 1044.72, '-137'],
[1232.17, 1007.96, '-123'],
[1037.54, 975.52, '-104'],
[1065.65, 703.04, '-107'],
[1361.93, 819.82, '-136'],
[1437.62, 592.75, '-144'],
[1225.68, 590.58, '-123'],
[1111.07, 417.58, '-111'],
[1080.79, 248.89, '-108'],
[1275.42, 264.03, '-128'],
[1420.32, 281.33, '-142'],
[1260.28, 456.5, '-126'],
[1344.63, 387.3, '-134'],
[909.95, 357.02, '-91'],  
]
# Convert the coordinates to a string format that can be inserted into the text widget
house_string = '\n'.join(str(coord) for coord in house_coordinates)
coordinates_text.insert(tk.END, house_string)


reproduce_button = tk.Button(root, text="Vis XY plassering", command=reproduce_drawing)
reproduce_button.grid(row=2, column=1, columnspan=1)


min_colormap_range_label = tk.Label(root, text="Min color map range:")
min_colormap_range_label.grid(row=5, column=0, sticky='e')
min_colormap_range_entry = tk.Entry(root, width=3)
min_colormap_range_entry.grid(row=5, column=1, sticky='w')

max_colormap_range_label = tk.Label(root, text="Max color map range:")
max_colormap_range_label.grid(row=6, column=0, sticky='e')
max_colormap_range_entry = tk.Entry(root, width=3)
max_colormap_range_entry.grid(row=6, column=1, sticky='w')

contour_lines_enabled = tk.BooleanVar(value=True)  # Checked by default
contour_lines_enabled_checkbutton = tk.Checkbutton(root, text="Enable contour lines", variable=contour_lines_enabled)
contour_lines_enabled_checkbutton.grid(row=9, column=0, columnspan=2, sticky='w') 

contour_lines_num_label = tk.Label(root, text="Antall linjer:")
contour_lines_num_label.grid(row=9, column=0, sticky='e')
contour_lines_num_entry = tk.Entry(root, width=3)
contour_lines_num_entry.grid(row=9, column=1, sticky='w') 
contour_lines_num_entry.insert(tk.END, "6")  # Default number of contour lines

contour_lines_width_label = tk.Label(root, text="Contour line width:")
contour_lines_width_label.grid(row=10, column=0, sticky='e')
contour_lines_width_entry = tk.Entry(root, width=3)
contour_lines_width_entry.grid(row=10, column=1, sticky='w')
contour_lines_width_entry.insert(tk.END, "1")  # Default contour line width


heatmap_button = tk.Button(root, text="Generate Heatmap", command=generate_heatmap)
heatmap_button.grid(row=20, column=0, columnspan=2, pady=15)

# Available colormaps
cmap_options = ["rosegold", "YlOrRd_r", "YlOrRd_r_opace", "viridis", "plasma", "inferno", "magma", "cividis", "brown_yellow_white", "brown_orange_white", 
                "brownred_orange_white", "brownred_yellow_white", "goldenrod", "rosegold chocolate", "rosegold chocolate_opace"]  # Add any other colormaps you want to include

# Create a StringVar to hold the current selection
selected_colormap = tk.StringVar(root)
selected_colormap.set(cmap_options[14])  # Set the default value

# Create the dropdown menu
cmap_menu = tk.OptionMenu(root, selected_colormap, *cmap_options)
cmap_menu.grid(row=19, column=0, columnspan=3)


root.mainloop()
