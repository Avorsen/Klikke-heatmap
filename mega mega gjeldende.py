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



def change_timeout():
    new_timeout = simpledialog.askinteger("Input", "Enter new timeout in milliseconds:", parent=root, minvalue=100, maxvalue=10000)
    if new_timeout is not None:  # Check if the user clicked OK and entered a value
        global custom_timeout
        custom_timeout = new_timeout
        average_z_display.set(f"Timeout set to {custom_timeout} ms")  # Optional: give feedback in the GUI


class TimedInputDialog(tk.Toplevel):
    def __init__(self, parent, title, prompt, timeout=1000):  # timeout now 1000 milliseconds
        super().__init__(parent)
        self.title(title)
        self.prompt = prompt
        self.timeout = timeout
        self.result = None
        self.timer_count = timeout  # Start counting from the full timeout in milliseconds
        self.create_widgets()
        self.update_timer()  # Start the countdown

    def create_widgets(self):
        self.label = tk.Label(self, text=self.prompt)
        self.label.pack(side="top", fill="x", padx=50, pady=10)

        self.timer_label = tk.Label(self, text=f"Auto-accepting in {self.timer_count} ms...")
        self.timer_label.pack(side="top", fill="x", padx=50, pady=5)
        
        self.entry = tk.Entry(self)
        self.entry.pack(padx=50, pady=10)
        self.entry.focus()
        
        self.button = tk.Button(self, text="OK", command=self.on_ok)
        self.button.pack(pady=10)

    def on_ok(self):
        self.result = self.entry.get()
        self.destroy()

    def update_timer(self):
        self.timer_count -= 100  # Update the countdown by reducing it by 100 ms
        if self.timer_count > 0:
            self.timer_label.config(text=f"Auto-accepting in {self.timer_count} ms...")
            self.after(100, self.update_timer)  # Schedule the next update every 100 ms
        else:
            self.on_ok()  # Time is up, close the dialog

def get_user_input(parent, title, prompt, timeout=1500):
    dialog = TimedInputDialog(parent, title, prompt, timeout)
    dialog.transient(parent)  # Make the dialog appear on top
    dialog.wait_window()  # Wait until the dialog is closed
    return dialog.result



class ClickHandler:
    def __init__(self, ax, text):
        self.ax = ax
        self.text = text
        self.cid = ax.figure.canvas.mpl_connect('button_press_event', self)
        self.points = []  # Store tuples of (plot object, coordinates, annotation)

    def __call__(self, event):
        if event.inaxes != self.ax:
            return

        x, y = round(event.xdata, 2), round(event.ydata, 2)

        # Handling different mouse buttons with Shift key consideration
        if event.button == 1:  # Left-click
            if event.key == 'shift':  # Check if Shift is pressed
                self.remove_point(x, y)
            else:
                z_value = get_user_input(self.ax.figure.canvas._tkcanvas, "Input", "Enter the Z value:", custom_timeout)  # 1000 ms timeout
                if z_value is not None and not z_value.startswith(('+', '-')):
                    z_value = '-' + z_value  # Ensure negative if not specified
                self.add_point(x, y, z_value)
        elif event.button == 3:  # Right-click
            z_value = get_user_input(self.ax.figure.canvas._tkcanvas, "Input", "Enter the Z value (positive):", custom_timeout)  # 1000 ms timeout
            if z_value is not None and not z_value.startswith(('+', '-')):
                z_value = '+' + z_value  # Automatically positive
            self.add_point(x, y, z_value)

    def add_point(self, x, y, z_value):
        if z_value is None:  # If user cancels the prompt
            return
        plot_obj = self.ax.plot(x, y, 'ko')[0]
        annotation = self.ax.annotate(z_value, (x, y), textcoords="offset points", xytext=(0,10), ha='center', color='red')
        self.points.append((plot_obj, (x, y, z_value), annotation))
        self.text.insert(tk.END, str((x, y, z_value)) + '\n')
        self.ax.figure.canvas.draw()

    def remove_point(self, x, y):
        # Find the closest point within some threshold
        closest_point = None
        min_dist = float('inf')
        for point in self.points:
            px, py, _ = point[1]
            dist = (px - x)**2 + (py - y)**2
            if dist < min_dist:
                min_dist = dist
                closest_point = point

        # Remove the point if it's close enough
        if closest_point and min_dist < 30:  # Threshold for "closeness"
            plot_obj, coords, annotation = closest_point
            plot_obj.remove()  # Remove the point from the plot
            annotation.remove()  # Remove the annotation
            self.points.remove(closest_point)  # Remove from the list
            self.update_text()  # Update the text widget
            self.ax.figure.canvas.draw()

    def update_text(self):
        # Refresh the text area
        self.text.delete('1.0', tk.END)
        for _, coords, _ in self.points:
            self.text.insert(tk.END, str(coords) + '\n')



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

    ax.invert_yaxis()  # Ensure the y-axis is inverted to match image orientation

    # Create a new top level window
    top = tk.Toplevel(root)
    top.title("Set Boundary")  # You can set title for the new window

    # Create the canvas on this new window
    canvas = FigureCanvasTkAgg(fig, master=top)
    canvas.draw()
    canvas.get_tk_widget().pack()  # Using pack() for simplicity


def upload_image():
    filepath = filedialog.askopenfilename()
    image = np.array(plt.imread(filepath))

    fig, ax = plt.subplots(figsize=(40,40))  # Define size of the figure here
    ax.imshow(image)


    global click_handler
    click_handler = ClickHandler(ax, coordinates_text)

    # Create a new top level window
    top = tk.Toplevel(root)
    top.title("Plantegning")  # You can set title for the new window

    # Create the canvas on this new window
    canvas = FigureCanvasTkAgg(fig, master=top)
    canvas.draw()
    canvas.get_tk_widget().pack()  # Using pack() for simplicity

    # Add updated instructions directly on the plot in Norwegian
    ax.text(0.01, 0.99, '1) Venstreklikk for å angi Z-verdi (blir negativ hvis ingen tegn angis).',
            transform=ax.transAxes, color='green', fontsize=14, va='top')
    ax.text(0.01, 0.97, '2) Høyreklikk for å angi Z-verdi (blir positiv hvis ingen tegn angis).',
            transform=ax.transAxes, color='green', fontsize=14, va='top')
    ax.text(0.01, 0.95, '3) Shift + Venstreklikk på et punkt for å fjerne det.',
            transform=ax.transAxes, color='green', fontsize=14, va='top')


def apply_z_adjustment():
    # Prompt user for Z-value adjustment
    adjustment = simpledialog.askfloat("Input", "Enter new Z-value adjustment:", parent=root)
    if adjustment is None:  # User cancelled the dialog
        return

    # Apply the adjustment to coordinates
    adjusted_coordinates = []
    raw_coordinates = coordinates_text.get("1.0", tk.END).strip().split('\n')
    
    # Process each coordinate line
    for coord in raw_coordinates:
        try:
            x, y, z = eval(coord)  # Safely evaluate the tuple
            adjusted_z = float(z) + adjustment  # Apply the adjustment
            adjusted_coordinates.append((x, y, adjusted_z))
        except (TypeError, ValueError, SyntaxError):
            average_z_display.set("Error in coordinate formatting. Ensure they are in '(x, y, z)' format.")
            return

    # Update the text box with adjusted coordinates
    coordinates_text.delete('1.0', tk.END)
    for coord in adjusted_coordinates:
        coordinates_text.insert(tk.END, str(coord) + '\n')
    average_z_display.set(f"Z-values adjusted by {adjustment}.")


# Create and register custom colormaps
colors = ["DarkOrange", "PaleGoldenrod", "white"]
cmap = mcolors.LinearSegmentedColormap.from_list("rosegold", colors)
plt.register_cmap(cmap=cmap)

colors2 = ["chocolate", "DarkOrange", "PaleGoldenrod", "white"]
cmap2 = mcolors.LinearSegmentedColormap.from_list("rosegold chocolate", colors2)
plt.register_cmap(cmap=cmap2)


colors_custom = ['#00204c', '#002080', '#008000']
cmap_custom = mcolors.LinearSegmentedColormap.from_list("custom_green_to_dark_blue", colors_custom)
plt.register_cmap(cmap=cmap_custom)

# Custom Blue Shades
colors_blue = ['#00204c', '#002080', 'lightblue']
cmap_custom_blue = mcolors.LinearSegmentedColormap.from_list("custom_blue", colors_blue)
plt.register_cmap(cmap=cmap_custom_blue)

# Assuming 'custom_blue' is already a registered colormap
original_cmap = plt.get_cmap("custom_blue")
colors_opaque = original_cmap(np.arange(original_cmap.N))

alpha = 0.84
colors_opaque[:, -1] = alpha 
new_cmap_name = "custom_blue_opaque"
new_cmap = mcolors.LinearSegmentedColormap.from_list(new_cmap_name, colors_opaque)
plt.register_cmap(cmap=new_cmap)



colors_custom = ['#00204c', '#002080', '#008000']
custom_cmap = mcolors.LinearSegmentedColormap.from_list("custom_green_to_dark_blue", colors_custom)
colors_array = custom_cmap(np.arange(custom_cmap.N))
alpha = 0.8  # Set opacity strength here
colors_array[:, -1] = alpha  # Change alpha values of all colors
opaque_cmap_custom = mcolors.LinearSegmentedColormap.from_list("custom_green_to_dark_blue_opaque", colors_array)
plt.register_cmap(cmap=opaque_cmap_custom)


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

def calculate_interpolation_data():
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

    # Use the selected interpolation method
    method = selected_interpolation_method.get()
    zi = scipy.interpolate.griddata((x, y), z, (xi, yi), method=method)

    boundary_path = None
    points2D = np.array(points)[:, :2]  # Only x and y for 2D hull

    if boundary_click_handler and boundary_click_handler.boundary_closed:
        boundary_path = Path(boundary_click_handler.boundary)
    else:
        # Create convex hull
        hull = ConvexHull(points2D)
        boundary_path = Path(hull.points[hull.vertices])

    # Create a mask for the region outside the contours
    points = np.column_stack((xi.flatten(), yi.flatten()))
    mask = ~boundary_path.contains_points(points).reshape(xi.shape)

    masked_zi = np.ma.array(zi, mask=mask)
    return x, y, z, xi, yi, masked_zi

def calculate_average_z():
    try:
        _, _, _, xi, yi, masked_zi = calculate_interpolation_data()

        # Calculate the area of each cell in the grid
        dx = np.diff(xi[0, :])[0]  # Difference between any two consecutive x values
        dy = np.diff(yi[:, 0])[0]  # Difference between any two consecutive y values
        cell_area = dx * dy

        # Calculate the total area that contributes to the weighted average
        total_area = np.sum(~masked_zi.mask) * cell_area

        # Calculate the weighted average Z
        total_weighted_z = np.sum(masked_zi * cell_area)
        weighted_average_z = total_weighted_z / total_area if total_area else np.nan  # Avoid division by zero

        average_z_display.set(f"Weighted Average Z: {weighted_average_z:.2f}")
    except Exception as e:
        average_z_display.set(f"Error: {str(e)}")

def visualize_weighted_calculation():
    try:
        _, _, _, xi, yi, masked_zi = calculate_interpolation_data()
        
        # Create a plot to show the interpolation and the cell areas
        fig, ax = plt.subplots(figsize=(10, 8))
        c = ax.pcolormesh(xi, yi, masked_zi, shading='auto', cmap='viridis')
        fig.colorbar(c, ax=ax, label='Interpolated Z values')
        
        # Optionally, highlight non-masked (used) areas
        ax.contourf(xi, yi, ~masked_zi.mask, 1, colors='none', hatches=['.'])

        ax.set_title("Interpolated Surface and Grid Cell Areas")
        ax.set_xlabel("X coordinate")
        ax.set_ylabel("Y coordinate")
        plt.show()
    except Exception as e:
        average_z_display.set(f"Error: {str(e)}")

def generate_heatmap():
    x, y, z, xi, yi, masked_zi = calculate_interpolation_data()

    # Ensure colormap range entries are set before using them
    min_colormap_range = float(min_colormap_range_entry.get()) if min_colormap_range_entry.get() else min(z) - 10
    max_colormap_range = float(max_colormap_range_entry.get()) if max_colormap_range_entry.get() else max(z)

    fig, ax = plt.subplots(figsize=(10,10))
    im = ax.imshow(masked_zi, extent=(min(x), max(x), min(y), max(y)), cmap=selected_colormap.get(), origin='lower', vmin=min_colormap_range, vmax=max_colormap_range)
    ax.invert_yaxis()  # Invert the y-axis for correct heatmap orientation

    # Optional background transparency
    fig.patch.set_visible(False)
    ax.axis('off')

    # Colorbar setup
    cbar = plt.colorbar(im, ax=ax, aspect=len(yi)//25, shrink=0.8)
    cbar.ax.set_ylabel('Setningsmålinger (mm)', rotation=-90, va="bottom")

    # Adding contour lines if enabled
    if contour_lines_enabled.get():
        contour_line_width = float(contour_lines_width_entry.get())
        num_contour_lines = int(contour_lines_num_entry.get())
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
    for x_coord, y_coord, z_value in zip(x, y, z):
        z_value_str = ('+' if z_value > 0 else '') + str(int(z_value)) if float(z_value).is_integer() else str(z_value)
        ax.text(x_coord, y_coord, z_value_str, fontsize=int(font_size.get()), ha='center', va='center', color=font_color.get())

    # Save the figure with z values
    plt.savefig(os.path.join(os.path.expanduser("~"), "Desktop", f"{saved_image_name}_with_z_values.png"), dpi=300, bbox_inches='tight', pad_inches=0)
    
    
    plt.show()
    print("Heatmap saved to Desktop. You can find the images with and without z-values, and dot on xy if ticked")




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
ascii_art_label.grid(row=0, column=1, columnspan=2)


upload_button = tk.Button(root, text="Velg Plantegning", command=upload_image)
upload_button.grid(row=1, column=0)

boundary_button = tk.Button(root, text="Set Boundary", command=set_boundary)
boundary_button.grid(row=2, column=0) 


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
contour_lines_num_entry.insert(tk.END, "10")  # Default number of contour lines

contour_lines_width_label = tk.Label(root, text="Contour line width:")
contour_lines_width_label.grid(row=10, column=0, sticky='e')
contour_lines_width_entry = tk.Entry(root, width=3)
contour_lines_width_entry.grid(row=10, column=1, sticky='w')
contour_lines_width_entry.insert(tk.END, "1")  # Default contour line width

font_size_label = tk.Label(root, text="Font Size:")
font_size_label.grid(row=11, column=0, sticky='e')
font_size = tk.StringVar(value="12")  # default font size
font_size_entry = tk.Entry(root, textvariable=font_size, width=3)
font_size_entry.grid(row=11, column=1, sticky='w')

# GUI setup for font color dropdown menu
color_options = ["black", "red", "blue", "green", "yellow", "purple"]
font_color_label = tk.Label(root, text="Font Color:")
font_color_label.grid(row=12, column=0, sticky='e')
font_color = tk.StringVar(value="black")  # default color (black)
color_menu = ttk.OptionMenu(root, font_color, font_color.get(), *color_options)
color_menu.grid(row=12, column=1, sticky='w')


cmap_label = tk.Label(root, text="Choose Colormap:")
cmap_label.grid(row=13, column=0, sticky='e')
cmap_options = ["rosegold chocolate", "rosegold chocolate_opace", "rosegold", "YlOrRd_r", "YlOrRd_r_opace", "viridis", "plasma", "inferno", "magma", "cividis",
                "custom_green_to_dark_blue", "custom_green_to_dark_blue_opaque", "custom_blue", "custom_blue_opaque"]

# Create a StringVar to hold the current selection
selected_colormap = tk.StringVar(root)
selected_colormap.set(cmap_options[1])  # Set the default value

# Create the colormap dropdown menu and position it
colormap_menu = ttk.OptionMenu(root, selected_colormap, selected_colormap.get(), *cmap_options)
colormap_menu.grid(row=13, column=1, sticky='w')

# List of interpolation methods
interpolation_methods = ['linear', 'nearest', 'cubic']
selected_interpolation_method = tk.StringVar(value=interpolation_methods[2])  # default to 'linear'

# Create the dropdown menu for interpolation methods
interpolation_method_label = tk.Label(root, text="Select Interpolation Method:")
interpolation_method_label.grid(row=14, column=0, sticky='e')
interpolation_method_menu = ttk.OptionMenu(root, selected_interpolation_method, selected_interpolation_method.get(), *interpolation_methods)
interpolation_method_menu.grid(row=14, column=1, sticky='w')

heatmap_button = tk.Button(root, text="Generate Heatmap", command=generate_heatmap)
heatmap_button.grid(row=22, column=0, columnspan=2, pady=15)


# Add an entry widget to display the average Z value
average_z_display = tk.StringVar(value="Average Z: Not calculated")
average_z_label = tk.Label(root, textvariable=average_z_display)
average_z_label.grid(row=20, column=1, columnspan=2, pady=5)

# Create a button to calculate the average Z value
calculate_z_button = tk.Button(root, text="Calculate Average Z", command=calculate_average_z)
calculate_z_button.grid(row=21, column=1, padx=5, pady=0)

# Add a button for visualizing the weighted Z calculation
visualize_button = tk.Button(root, text="Visualize Weighted Calculation", command=visualize_weighted_calculation)
visualize_button.grid(row=21, column=0, padx=5, pady=0)

# UI components for Z adjustment
# Button to trigger Z adjustment dialog
adjustment_button = tk.Button(root, text="Adjust Z Values", command=apply_z_adjustment)
adjustment_button.grid(row=2, column=1)


timeout_button = tk.Button(root, text="Change Timeout", command=change_timeout)
timeout_button.grid(row=1, column=1)  # Adjust positioning as needed

custom_timeout = 1000  # Default timeout value

root.mainloop()
