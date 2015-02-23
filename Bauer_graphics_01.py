#Bauer, Matt
#1000-631-613
#2015-02-8
#Assignment_03

import math

class ClWorld:
    def __init__(self, objects=[], canvases=[]):
        self.vertices = []
        self.faces = []
        self.vrp = None
        self.vpn = None
        self.vup = None
        self.prp = None
        self.window = None
        self.viewport = None
        self.center_of_window = None
        self.composite_transform = [[1, 0, 0, 0],
                                    [0, 1, 0, 0],
                                    [0, 0, 1, 0],
                                    [0, 0, 0, 1]]
        self.object_origin = [0.0, 0.0, 0.0, 1.0]

        # mouse variables
        self._x = None
        self._y = None

        self.objects = objects
        self.canvases = canvases

    def add_canvas(self, canvas):
        self.canvases.append(canvas)
        canvas.world = self

    def create_graphic_objects(self, canvas, event=None):

        # track the origin of the object being displayed
        self.object_origin = multiply_vector(self.composite_transform, self.object_origin)
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


        scale_t = scale_matrix(max(self.window[0], self.window[1])-min(self.window[0], self.window[1]),
                               max(self.window[2], self.window[3])-min(self.window[2], self.window[3]),
                               max(self.window[4], self.window[5])-min(self.window[4], self.window[5]))
        t2_1 = translation_matrix(-min(self.window[0], self.window[1]),
                                  -min(self.window[2], self.window[3]),
                                  -min(self.window[4], self.window[5]))
        window_transform = [[1, 0, 0, -self.window[0]],
                            [0, -1, 0, self.window[3]],
                            [0, 0, 1, 0],
                            [0, 0, 0, 1]]
        scale_transform = scale_matrix(1/size_x * ratio_x, 1/size_y * ratio_y, 1)
        viewport_transform = translation_matrix(-tran_x, -tran_y, 0)

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

    def move(self, move_amount):
        t1 = translation_matrix(-move_amount[0], -move_amount[1], -move_amount[2])
        self.composite_transform = multiply_matrix(t1, self.composite_transform)

    def unit_cube_clip(self):
        t1 = translation_matrix(self.vrp[0], self.vrp[1], self.vrp[2])

        rx = rotate_to_plane_matrix('x', self.vpn[0], self.vpn[1], self.vpn[2])
        vpn_prime = multiply_vector(rx, self.vpn)

        ry = rotate_to_plane_matrix('y', vpn_prime[0], vpn_prime[1], vpn_prime[2])

        vup_prime = multiply_vector(rx, self.vup)
        vup_prime = multiply_vector(ry, vup_prime)
        rz = rotate_to_plane_matrix('z', vup_prime[0], vup_prime[1], vup_prime[2])

        direction_of_projection = [self.center_of_window[0] - self.prp[0],
                                   self.center_of_window[1] - self.prp[1],
                                   -self.prp[2]]

        shear = shear_matrix(direction_of_projection[0], direction_of_projection[1], direction_of_projection[2])

        t2 = translation_matrix(min(self.window[0], self.window[1]),
                                min(self.window[2], self.window[3]),
                                min(self.window[4], self.window[5]))

        scale = scale_matrix(1/(max(self.window[0], self.window[1])-min(self.window[0], self.window[1])),
                             1/(max(self.window[2], self.window[3])-min(self.window[2], self.window[3])),
                             1/(max(self.window[4], self.window[5])-min(self.window[4], self.window[5])))

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

    def scale(self, s_x, s_y, s_z, point=None):
        if point:
            t1 = translation_matrix(point[0], point[1], point[2])
            t2 = translation_matrix(-point[0], -point[1], -point[2])
        else:
            t1 = [[1, 0, 0, 0],
                  [0, 1, 0, 0],
                  [0, 0, 1, 0],
                  [0, 0, 0, 1]]
            t2 = [[1, 0, 0, 0],
                  [0, 1, 0, 0],
                  [0, 0, 1, 0],
                  [0, 0, 0, 1]]

        scale_matrix = [[s_x, 0, 0, 0],
                        [0, s_y, 0, 0],
                        [0, 0, s_z, 0],
                        [0, 0, 0, 1]]
        self.composite_transform = multiply_matrix(t2, multiply_matrix(scale_matrix, multiply_matrix(t1, self.composite_transform)))

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

'''
************************************************************

************************************************************
'''


def clip(point_a, point_b):
    a_byte_array = make_byte_array(point_a)
    b_byte_array = make_byte_array(point_b)

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


def scale_matrix(s_x, s_y, s_z):
    return [[s_x, 0, 0, 0],
            [0, s_y, 0, 0],
            [0, 0, s_z, 0],
            [0, 0, 0, 1]]


def shear_matrix(a, b, c):
    if c == 0:
        return [[1, 0, 0, 0],
                [0, 1, 0, 0],
                [0, 0, 1, 0],
                [0, 0, 0, 1]]

    return [[1, 0, -a/c, 0],
            [0, 1, -b/c, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1]]


def translation_matrix(t_x, t_y, t_z):
    return [[1, 0, 0, -t_x],
            [0, 1, 0, -t_y],
            [0, 0, 1, -t_z],
            [0, 0, 0, 1]]


def rotate_to_plane_matrix(axis_of_rotation, a, b, c):
    if axis_of_rotation == 'x':
        denominator = math.sqrt(b*b+c*c)
        if denominator == 0:
            return [[1, 0, 0, 0],
                    [0, 1, 0, 0],
                    [0, 0, 1, 0],
                    [0, 0, 0, 1]]
        cos = c / denominator
        sin = b / denominator
        return [[1, 0, 0, 0],
                [0, cos, -sin, 0],
                [0, sin, cos, 0],
                [0, 0, 0, 1]]
    elif axis_of_rotation == 'y':
        denominator = math.sqrt(a*a+c*c)
        if denominator == 0:
            return [[1, 0, 0, 0],
                    [0, 1, 0, 0],
                    [0, 0, 1, 0],
                    [0, 0, 0, 1]]
        cos = c / denominator
        sin = a / denominator
        return [[cos, 0, -sin, 0],
                [0, 1, 0, 0],
                [sin, 0, cos, 0],
                [0, 0, 0, 1]]
    elif axis_of_rotation == 'z':
        denominator = math.sqrt(a*a+b*b)
        if denominator == 0:
            return [[1, 0, 0, 0],
                    [0, 1, 0, 0],
                    [0, 0, 1, 0],
                    [0, 0, 0, 1]]
        cos = b / denominator
        sin = a / denominator
        return [[cos, -sin, 0, 0],
                [sin, cos, 0, 0],
                [0, 0, 1, 0],
                [0, 0, 0, 1]]

    return [[1, 0, 0, 0],
            [0, 1, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1]]
