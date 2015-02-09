#Bauer, Matt
#1000-631-613
#2015-02-8
#Assignment_01

import math

class ClWorld:
    def __init__(self, objects=[], canvases=[]):
        self.vertices = []
        self.faces = []
        self.window = None
        self.viewport = None

        # mouse variables
        self._x = None
        self._y = None

        self.objects = objects
        self.canvases = canvases

    def add_canvas(self, canvas):
        self.canvases.append(canvas)
        canvas.world = self

    def create_graphic_objects(self, canvas, event=None):
        if event:
            if self._y:
                y_rotation = 0
                if self._y > event.y:
                    y_rotation = self._y - event.y
                elif self._y < event.y:
                    y_rotation = self._y - event.y

                x_rotation = 0
                if self._x > event.x:
                    x_rotation = self._x - event.x
                elif self._x < event.x:
                    x_rotation = self._x - event.x

                composite_matrix = None
                y_rotation_matrix = None
                if y_rotation != 0:
                    # rotate around z axis
                    y_rotation_matrix = [[1, 0, 0, 0],
                                         [0, math.cos(y_rotation*.001), -math.sin(y_rotation*.001), 0],
                                         [0, math.sin(y_rotation*.001), math.cos(y_rotation*.001), 0],
                                         [0, 0, 0, 1]]
                    composite_matrix = y_rotation_matrix

                if x_rotation != 0:
                    # rotate around y axis
                    x_rotation_matrix = [[math.cos(x_rotation*.001), 0, math.sin(x_rotation*.001), 0],
                                         [0, 1, 0, 0],
                                         [-math.sin(x_rotation*.001), 0, math.cos(x_rotation*.001), 0],
                                         [0, 0, 0, 1]]

                    if y_rotation_matrix:
                        composite_matrix = multiply_matrix(y_rotation_matrix, x_rotation_matrix)
                    else:
                        composite_matrix = x_rotation_matrix

                if composite_matrix:
                    for x in range(0, len(self.vertices)):
                        self.vertices[x] = multiply_vector(composite_matrix, self.vertices[x])

            # set local x and y values for the next method run
            self._x = event.x
            self._y = event.y
        else:
            self._y = None

        min_x = float(canvas.cget("width")) * self.viewport[0]
        max_x = float(canvas.cget("width")) * self.viewport[2]

        min_y = float(canvas.cget("height")) * self.viewport[1]
        max_y = float(canvas.cget("height")) * self.viewport[3]

        #print("screen min max")
        #print(min_x)
        #print(max_x)
        #print(min_y)
        #print(max_y)

        self.objects.append(canvas.create_rectangle(min_x, max_y, max_x, min_y))

        ratio_x = (max_x - min_x)
        tran_x = min_x

        ratio_y = (max_y - min_y)
        tran_y = min_y

        size_x = self.window[2] - self.window[0]
        size_y = self.window[3] - self.window[1]

        window_transform = [[1, 0, 0, -self.window[0]],
                            [0, -1, 0, self.window[3]],
                            [0, 0, 1, 0],
                            [0, 0, 0, 1]]

        scale_transform = [[1/size_x * ratio_x, 0, 0, 0],
                           [0, 1/size_y * ratio_y, 0, 0],
                           [0, 0, 1, 0],
                           [0, 0, 0, 1]]

        viewport_transform = [[1, 0, 0, tran_x],
                              [0, 1, 0, tran_y],
                              [0, 0, 1, 0],
                              [0, 0, 0, 1]]

        composite_transform = multiply_matrix(viewport_transform, multiply_matrix(scale_transform, window_transform))

        new_verts = self.vertices.copy()
        for x in range(0, len(new_verts)):
            new_verts[x] = multiply_vector(composite_transform, new_verts[x])
            pass

        for face in self.faces:
            points = []
            for x in range(0, len(face)):
                points.append(new_verts[face[x]][0])
                points.append(new_verts[face[x]][1])

            self.objects.append(canvas.create_polygon(points, fill="red"))
            self.objects.append(canvas.create_line(points))

    def redisplay(self, canvas, event=None):
        canvas.delete("all")
        self.create_graphic_objects(canvas, event)

    def reset_lists(self):
        self.vertices = []
        self.faces = []
        self.window = None
        self.viewport = None


def multiply_vector(matrix, vector):
    new_vector = []
    for y in range(0, 4):
        new_value = 0
        for x in range(0, 4):

            new_value += matrix[y][x]*vector[x]
        new_vector.append(new_value)

    return new_vector


def multiply_matrix(matrix_a, matrix_b):
    result = [[0, 0, 0, 0],
              [0, 0, 0, 0],
              [0, 0, 0, 0],
              [0, 0, 0, 0]]
    # iterate through rows of X
    for i in range(len(matrix_a)):
        # iterate through columns of Y
        for j in range(len(matrix_b[0])):
            # iterate through rows of Y
            for k in range(len(matrix_b)):
                result[i][j] += matrix_a[i][k] * matrix_b[k][j]
    return result
