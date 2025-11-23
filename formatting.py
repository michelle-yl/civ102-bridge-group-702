import mplcursors


# shows coordinate values when hovering over plot with cursor
def format_cursor(points):
    # Parameters:
        # points (list): all coordinate points of plot
    # Returns:
        # cursor (interactive): shows data points when used with matplotlib plot


    # show cursor on hover
    cursor = mplcursors.cursor(points, hover=True)


    @cursor.connect("add")
    def on_add(sel):
        x, y = sel.target
       
        # format cursor coordinate display
        sel.annotation.set(text=f"x = {x:.0f}\ny = {y:.4f}")
        sel.annotation.arrow_patch.set(visible=False)
        sel.annotation.get_bbox_patch().set(color="lightgrey")


    return cursor