import mplcursors
from matplotlib.patches import ConnectionStyle

def format_cursor(points):
    cursor = mplcursors.cursor(points, hover=True)

    @cursor.connect("add")
    def on_add(sel):
        x, y = sel.target
        
        sel.annotation.set(text=f"x = {x:.0f}\ny = {y:.2f}")
        sel.annotation.arrow_patch.set(visible=False)
        sel.annotation.get_bbox_patch().set(color="lightgrey")

    return cursor