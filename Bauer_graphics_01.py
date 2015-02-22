#Bauer, Matt
#1000-631-613
#2015-02-8
#Assignment_03

import math

class ClWorld:
    def __init__(self, objects=[], canvases=[]):
        self.vertices = []
        self.faces = []
        self.r = None
        self.n = None
        self.u = None
        self.p = None
        self.window = None
        self.viewport = None
        self.center_of_window = None
        self.composite_transform = [[1, 0, 0, 0],
                                    [0, 1, 0, 0],
                                    [0, 0, 1, 0],
                                    [0, 0, 0, 1]]

        # mouse variables
        self._x = None
        self._y = None

        self.objects = objects
        self.canvases = canvases

    def add_canvas(self, canvas):
        self.canvases.append(canvas)
        canvas.world = self

    def create_graphic_objects(self, canvas, event=None):
        for x in range(0, len(self.vertices)):
            self.vertices[x] = multiply_vector(self.composite_transform, self.vertices[x])

        min_x = float(canvas.cget("width")) * self.viewport[0]
        max_x = float(canvas.cget("width")) * self.viewport[2]

        min_y = float(canvas.cget("height")) * self.viewport[1]
        max_y = float(canvas.cget("height")) * self.viewport[3]

        self.objects.append(canvas.create_rectangle(min_x, max_y, max_x, min_y))

        ratio_x = (max_x - min_x)
        tran_x = min_x

        ratio_y = (max_y - min_y)
        tran_y = min_y

        size_x = self.window[1] - self.window[0]
        size_y = self.window[3] - self.window[2]

        scale_t = [[(max(self.window[0], self.window[1])-min(self.window[0], self.window[1])), 0, 0, 0],
                    [0, (max(self.window[2], self.window[3])-min(self.window[2], self.window[3])), 0, 0],
                    [0, 0, (max(self.window[4], self.window[5])-min(self.window[4], self.window[5])), 0],
                    [0, 0, 0, 1]]

        t2_1 = [[1, 0, 0, min(self.window[0], self.window[1])],
               [0, 1, 0, min(self.window[2], self.window[3])],
               [0, 0, 1, min(self.window[4], self.window[5])],
               [0, 0, 0, 1]]

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

        composite_transform = multiply_matrix(viewport_transform, multiply_matrix(scale_transform, multiply_matrix(window_transform, multiply_matrix(t2_1, scale_t))))

        lines = self.unit_cube_clip()
        # print(lines)
        for line_index in range(0, len(lines)):
            # print(lines[line_index])
            for point_index in range(0, len(lines[line_index])):
                lines[line_index][point_index].append(1.0)
                lines[line_index][point_index] = multiply_vector(composite_transform, lines[line_index][point_index])

            self.objects.append(canvas.create_line(lines[line_index][0][0], lines[line_index][0][1],
                                                   lines[line_index][1][0], lines[line_index][1][1]))

        self.composite_transform = [[1, 0, 0, 0],
                                    [0, 1, 0, 0],
                                    [0, 0, 1, 0],
                                    [0, 0, 0, 1]]

    def unit_cube_clip(self):
        t1 = [[1, 0, 0, -self.r[0]],
              [0, 1, 0, -self.r[1]],
              [0, 0, 1, -self.r[2]],
              [0, 0, 0, 1]]

        denom_rx = math.sqrt(self.n[1] * self.n[1] + self.n[2] * self.n[2])
        if denom_rx == 0:
            rx = [[1, 0, 0, 0],
                  [0, 1, 0, 0],
                  [0, 0, 1, 0],
                  [0, 0, 0, 1]]
        else:
            rx = [[1, 0, 0, 0],
                  [0, self.n[2]/denom_rx,  -self.n[1]/denom_rx, 0],
                  [0, self.n[1]/denom_rx,  self.n[2]/denom_rx, 0],
                  [0, 0, 0, 1]]

        denom_ry = math.sqrt(self.n[0] * self.n[0] + self.n[2] * self.n[2] + self.n[1]*self.n[1])
        if denom_ry == 0:
            ry = [[1, 0, 0, 0],
                  [0, 1, 0, 0],
                  [0, 0, 1, 0],
                  [0, 0, 0, 1]]
        else:
            ry = [[denom_rx/denom_ry, 0, -self.n[0]/denom_ry, 0],
                  [0, 1, 0, 0],
                  [self.n[0]/denom_ry, 0, denom_rx/denom_ry, 0],
                  [0, 0, 0, 1]]

        vup = self.u.copy()
        vup.append(1)
        vup = multiply_vector(rx, vup)
        vup = multiply_vector(ry, vup)

        denom_rz = math.sqrt(vup[0]*vup[0] + vup[1]*vup[1])
        if denom_rz == 0:
            rz = [[1, 0, 0, 0],
                  [0, 1, 0, 0],
                  [0, 0, 1, 0],
                  [0, 0, 0, 1]]
        else:
            rz = [[vup[1]/denom_rz, -vup[0]/denom_rz, 0, 0],
                  [vup[0]/denom_rz, vup[1]/denom_rz, 0, 0],
                  [0, 0, 1, 0],
                  [0, 0, 0, 1]]

        direction_of_projection = [self.center_of_window[0] - self.p[0],
                                   self.center_of_window[1] - self.p[1],
                                   - self.p[2]]

        shear = [[1, 0, -direction_of_projection[0]/direction_of_projection[2], 0],
                 [0, 1, -direction_of_projection[1]/direction_of_projection[2], 0],
                 [0, 0, 1, 0],
                 [0, 0, 0, 1]]

        t2 = [[1, 0, 0, -min(self.window[0], self.window[1])],
              [0, 1, 0, -min(self.window[2], self.window[3])],
              [0, 0, 1, -min(self.window[4], self.window[5])],
              [0, 0, 0, 1]]

        scale = [[1/(max(self.window[0], self.window[1])-min(self.window[0], self.window[1])), 0, 0, 0],
                 [0, 1/(max(self.window[2], self.window[3])-min(self.window[2], self.window[3])), 0, 0],
                 [0, 0, 1/(max(self.window[4], self.window[5])-min(self.window[4], self.window[5])), 0],
                 [0, 0, 0, 1]]

        composite = multiply_matrix(scale, multiply_matrix(t2, multiply_matrix(shear, multiply_matrix(rz, multiply_matrix(ry , multiply_matrix(rx, t1))))))

        verts = self.vertices.copy()
        for x in range(0, len(verts)):
            verts[x] = multiply_vector(composite, verts[x])

        lines = []
        for face in self.faces:
            points = []
            for x in range(0, len(face)):
                points.append([verts[face[x]][0],
                               verts[face[x]][1],
                               verts[face[x]][2],
                               1.0])

            for x in range(0, len(points)):
                if x == len(points)-1:
                    line = clip(points[x], points[0])
                else:
                    line = clip(points[x], points[x + 1])

                if line:
                    lines.append(line)
        return lines

    def rotate_around_a_line(self, point_a, point_b, theta):
        point_a[0] = float(point_a[0])
        point_a[1] = float(point_a[1])
        point_a[2] = float(point_a[2])

        point_b[0] = float(point_b[0])
        point_b[1] = float(point_b[1])
        point_b[2] = float(point_b[2])
        point_b.append(1.0)

        theta = math.radians(theta)

        transpose = [[1, 0, 0, -point_a[0]],
                     [0, 1, 0, -point_a[1]],
                     [0, 0, 1, -point_a[2]],
                     [0, 0, 0, 1]]
        point_b = multiply_vector(transpose, point_b)
        # NOTE: use point_b' for other matrices and continue using the updated values

        r_x_denominator = math.sqrt(point_b[2]*point_b[2] + point_b[1]*point_b[1])
        r_x = [[1, 0, 0, 0],
               [0, point_b[2]/r_x_denominator, -point_b[1]/r_x_denominator, 0],
               [0, point_b[1]/r_x_denominator, point_b[2]/r_x_denominator, 0],
               [0, 0, 0, 1]]

        r_y_denominator = math.sqrt(point_b[2]*point_b[2] + point_b[0]*point_b[0] + point_b[1]*point_b[1])
        r_y = [[r_x_denominator/r_y_denominator, 0, -point_b[0]/r_y_denominator, 0],
               [0, 1, 0, 0],
               [point_b[0]/r_y_denominator, 0, r_x_denominator/r_y_denominator, 0],
               [0, 0, 0, 1]]

        rotation = [[math.cos(theta), -math.sin(theta), 0, 0],
                    [math.sin(theta), math.cos(theta), 0, 0],
                    [0, 0, 1, 0],
                    [0, 0, 0, 1]]

        r_y_transpose = [[r_x_denominator/r_y_denominator, 0, point_b[0]/r_y_denominator, 0],
                         [0, 1, 0, 0],
                         [-point_b[0]/r_y_denominator, 0, r_x_denominator/r_y_denominator, 0],
                         [0, 0, 0, 1]]

        r_x_transpose = [[1, 0, 0, 0],
                         [0, point_b[2]/r_x_denominator, point_b[1]/r_x_denominator, 0],
                         [0, -point_b[1]/r_x_denominator, point_b[2]/r_x_denominator, 0],
                         [0, 0, 0, 1]]

        transpose_inverse = [[1, 0, 0, point_a[0]],
                             [0, 1, 0, point_a[1]],
                             [0, 0, 1, point_a[2]],
                             [0, 0, 0, 1]]

        self.composite_transform = multiply_matrix(rotation, multiply_matrix(
            r_y, multiply_matrix(r_x, multiply_matrix(transpose, self.composite_transform))))
        self.composite_transform = multiply_matrix(transpose_inverse, multiply_matrix(
            r_x_transpose, multiply_matrix(r_y_transpose, self.composite_transform)))

    def rotate_theta(self, axis, theta):
        theta = math.radians(theta)
        rotation = [[1, 0, 0, 0],
                    [0, 1, 0, 0],
                    [0, 0, 1, 0],
                    [0, 0, 0, 1]]
        if axis == 1:
            rotation = [[1, 0, 0, 0],
                        [0, math.cos(theta), -math.sin(theta), 0],
                        [0, math.sin(theta), math.cos(theta), 0],
                        [0, 0, 0, 1]]
        elif axis == 2:
            rotation = [[math.cos(theta), 0, math.sin(theta), 0],
                        [0, 1, 0, 0],
                        [-math.sin(theta), 0, math.cos(theta), 0],
                        [0, 0, 0, 1]]
        elif axis == 3:
            rotation = [[math.cos(theta), -math.sin(theta), 0, 0],
                        [math.sin(theta), math.cos(theta), 0, 0],
                        [0, 0, 1, 0],
                        [0, 0, 0, 1]]
        self.composite_transform = multiply_matrix(rotation, self.composite_transform)

    def scale(self, s_x, s_y, s_z):
        scale_matrix = [[s_x, 0, 0, 0],
                        [0, s_y, 0, 0],
                        [0, 0, s_z, 0],
                        [0, 0, 0, 1]]
        self.composite_transform = multiply_matrix(scale_matrix, self.composite_transform)

    def redisplay(self, canvas, event=None):
        canvas.delete("all")
        self.create_graphic_objects(canvas, event)

    def reset_lists(self):
        self.vertices = []
        self.faces = []
        self.window = None
        self.viewport = None

    def rotate_with_mouse(self, event):
        if event:
            if self._x:
                y_rotation = 0
                if self._y > event.y:
                    y_rotation = -1
                elif self._y < event.y:
                    y_rotation = 1

                x_rotation = 0
                if self._x > event.x:
                    x_rotation = -1
                elif self._x < event.x:
                    x_rotation = 1

                if y_rotation != 0:
                    self.rotate_theta(1, 10 * y_rotation)

                if x_rotation != 0:
                    self.rotate_theta(2, 10 * x_rotation)

            # set local x and y values for the next method run
            self._x = event.x
            self._y = event.y
        else:
            self._x = None
            self._y = None

    def scale_with_mouse(self, event):
        if event:
            if self._y:
                y_rotation = 0
                if self._y > event.y:
                    y_rotation = -1
                elif self._y < event.y:
                    y_rotation = 1

                if y_rotation != 0:
                    self.scale(1.0 + y_rotation*.1, 1.0 + y_rotation*.1, 1.0 + y_rotation*.1)

            # set local x and y values for the next method run
            self._y = event.y
        else:
            self._y = None


def clip(point_a, point_b):
    a_byte_array = make_byte_array(point_a)
    b_byte_array = make_byte_array(point_b)

    and_byte_array = []
    or_byte_array = []

    or_bytes = False

    for x in range(0, len(a_byte_array)):
        if a_byte_array[x] and b_byte_array[x]:
            return None

        if a_byte_array[x] or b_byte_array[x]:
            or_bytes = True

    if not or_bytes:
        return [point_a, point_b]

    new_points = []
    move_vector = subtract_vectors(point_b, point_a)

    plane_coefficients = [[1, 0, 0, -1],
                          [1, 0, 0, 0],
                          [0, 1, 0, -1],
                          [0, 1, 0, 0],
                          [0, 0, 1, -1],
                          [0, 0, 1, 0]]

    if 0 <= point_a[0] <= 1 and 0 <= point_a[1] <= 1 and 0 <= point_a[2] <= 1:
        new_points.append(point_a)
    if 0 <= point_b[0] <= 1 and 0 <= point_b[1] <= 1 and 0 <= point_b[2] <= 1:
        new_points.append(point_b)

    for plane_index in range(0, len(plane_coefficients)):
        t = -(plane_coefficients[plane_index][0] * point_a[0] +
              plane_coefficients[plane_index][1] * point_a[1] +
              plane_coefficients[plane_index][2] * point_a[2] +
              plane_coefficients[plane_index][3]) / \
             (plane_coefficients[plane_index][0] * move_vector[0] +
              plane_coefficients[plane_index][1] * move_vector[1] +
              plane_coefficients[plane_index][2] * move_vector[2] + .000000000000000000000000001)
        if 0 <= t <= 1:
            # calculate new x y z point
            point = [move_vector[0] * t + point_a[0],
                     move_vector[1] * t + point_a[1],
                     move_vector[2] * t + point_a[2],
                     1.0]
            if 0 <= point[0] <= 1 and 0 <= point[1] <= 1 and 0 <= point[2] <= 1:
                new_points.append(point)
                if len(new_points) == 2:
                    return new_points
    return None


def transpose(matrix4):
    transpose_mat = [[matrix4[0][0], matrix4[1][0], matrix4[2][0], matrix4[3][0]],
                     [matrix4[0][1], matrix4[1][1], matrix4[2][1], matrix4[3][1]],
                     [matrix4[0][2], matrix4[1][2], matrix4[2][2], matrix4[3][2]],
                     [matrix4[0][3], matrix4[1][3], matrix4[2][3], matrix4[3][3]]]
    return transpose_mat


def make_byte_array(vector):

    byte_array = []

    for x in range(0, len(vector)):
        if vector[x] > 1:
            byte_array.append(True)
        else:
            byte_array.append(False)

        if vector[x] < 0:
            byte_array.append(True)
        else:
            byte_array.append(False)
    return byte_array


def subtract_vectors(vec_a, vec_b):
    return [vec_a[0] - vec_b[0],
            vec_a[1] - vec_b[1],
            vec_a[2] - vec_b[2]]


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
