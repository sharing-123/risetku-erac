import struct
from pathlib import Path
import matplotlib
import numpy as np
from PIL import Image, ImageDraw, ImageFont

ASSETS = Path(__file__).resolve().parent.parent / "assets"
ASSETS.mkdir(exist_ok=True)
STYLE = {"C": (0.68, (92, 108, 128)), "N": (0.66, (52, 110, 235)), "O": (0.64, (232, 76, 61)), "H": (0.40, (232, 236, 242))}


def polar(origin, length, degrees, z=0.0):
    angle = np.radians(degrees)
    return np.array([origin[0] + length * np.cos(angle), origin[1] + length * np.sin(angle), z])


def paracetamol():
    centre = np.zeros(3)
    atoms = [("C", polar(centre, 1.39, a)) for a in range(0, 360, 60)]
    bonds = [(i, (i + 1) % 6) for i in range(6)]

    def attach(parent, element, position):
        atoms.append((element, position))
        bonds.append((parent, len(atoms) - 1))
        return len(atoms) - 1

    hydroxyl = attach(3, "O", polar(atoms[3][1], 1.36, 180))
    attach(hydroxyl, "H", polar(atoms[hydroxyl][1], 0.96, 109))
    nitrogen = attach(0, "N", polar(atoms[0][1], 1.41, 0))
    attach(nitrogen, "H", polar(atoms[nitrogen][1], 1.01, -60))
    carbonyl = attach(nitrogen, "C", polar(atoms[nitrogen][1], 1.36, 55))
    attach(carbonyl, "O", polar(atoms[carbonyl][1], 1.23, 115))
    methyl = attach(carbonyl, "C", polar(atoms[carbonyl][1], 1.50, -5))
    attach(methyl, "H", polar(atoms[methyl][1], 1.09, 65))
    for z in (0.89, -0.89):
        attach(methyl, "H", polar(atoms[methyl][1], 0.63, -55, z))
    for ring, angle in [(1, 60), (2, 120), (4, 240), (5, 300)]:
        attach(ring, "H", polar(atoms[ring][1], 1.08, angle))
    elements = [a[0] for a in atoms]
    xyz = np.array([a[1] for a in atoms])
    return elements, xyz - xyz.mean(0), bonds


def icosphere(level=3):
    t = (1 + 5 ** 0.5) / 2
    vertices = [np.array(v, float) / np.linalg.norm(v) for v in [(-1, t, 0), (1, t, 0), (-1, -t, 0), (1, -t, 0), (0, -1, t), (0, 1, t),
                (0, -1, -t), (0, 1, -t), (t, 0, -1), (t, 0, 1), (-t, 0, -1), (-t, 0, 1)]]
    faces = [(0, 11, 5), (0, 5, 1), (0, 1, 7), (0, 7, 10), (0, 10, 11), (1, 5, 9), (5, 11, 4), (11, 10, 2), (10, 7, 6), (7, 1, 8),
             (3, 9, 4), (3, 4, 2), (3, 2, 6), (3, 6, 8), (3, 8, 9), (4, 9, 5), (2, 4, 11), (6, 2, 10), (8, 6, 7), (9, 8, 1)]
    triangles = [[vertices[i] for i in f] for f in faces]
    unit = lambda v: v / np.linalg.norm(v)
    for _ in range(level):
        finer = []
        for a, b, c in triangles:
            ab, bc, ca = unit(a + b), unit(b + c), unit(c + a)
            finer += [[a, ab, ca], [b, bc, ab], [c, ca, bc], [ab, bc, ca]]
        triangles = finer
    return np.array(triangles)


def cylinder(start, end, radius, sides=24):
    axis = (end - start) / np.linalg.norm(end - start)
    helper = np.array([0.0, 0.0, 1.0]) if abs(axis[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    u = np.cross(axis, helper)
    u /= np.linalg.norm(u)
    v = np.cross(axis, u)
    ring = [radius * (np.cos(a) * u + np.sin(a) * v) for a in np.linspace(0, 2 * np.pi, sides + 1)]
    triangles = []
    for p, q in zip(ring[:-1], ring[1:]):
        triangles += [[start + p, start + q, end + p], [start + q, end + q, end + p]]
    return np.array(triangles)


def write_stl(path, elements, xyz, bonds):
    sphere = icosphere()
    mesh = [sphere * STYLE[e][0] * 0.62 + p for e, p in zip(elements, xyz)]
    mesh += [cylinder(xyz[i], xyz[j], 0.13) for i, j in bonds]
    triangles = np.concatenate(mesh)
    with open(path, "wb") as handle:
        handle.write(b"paracetamol".ljust(80, b" ") + struct.pack("<I", len(triangles)))
        for t in triangles:
            normal = np.cross(t[1] - t[0], t[2] - t[0])
            normal = normal / (np.linalg.norm(normal) or 1)
            handle.write(struct.pack("<12fH", *normal, *t.ravel(), 0))


def shaded_sphere(radius, colour):
    y, x = np.mgrid[-radius:radius + 1, -radius:radius + 1] / radius
    r2 = x ** 2 + y ** 2
    z = np.sqrt(np.clip(1 - r2, 0, 1))
    light = np.array([-0.5, -0.6, 0.62])
    light /= np.linalg.norm(light)
    diffuse = np.clip(x * light[0] + y * light[1] + z * light[2], 0, 1)
    specular = diffuse ** 28
    rgb = np.array(colour)[None, None, :] * (0.22 + 0.82 * diffuse[..., None]) + 235 * specular[..., None]
    alpha = np.clip((1 - np.sqrt(r2)) * radius, 0, 1) * 255
    return Image.fromarray(np.dstack([np.clip(rgb, 0, 255), alpha]).astype(np.uint8), "RGBA")


def rotation(ax, ay, az):
    a, b, c = np.radians([ax, ay, az])
    rx = np.array([[1, 0, 0], [0, np.cos(a), -np.sin(a)], [0, np.sin(a), np.cos(a)]])
    ry = np.array([[np.cos(b), 0, np.sin(b)], [0, 1, 0], [-np.sin(b), 0, np.cos(b)]])
    rz = np.array([[np.cos(c), -np.sin(c), 0], [np.sin(c), np.cos(c), 0], [0, 0, 1]])
    return rz @ ry @ rx


def banner(path, elements, xyz, bonds, width=1920, height=640, scale=2):
    w, h = width * scale, height * scale
    gx, gy = np.meshgrid(np.linspace(0, 1, w), np.linspace(0, 1, h))
    glow = np.exp(-(((gx - 0.74) * 2.4) ** 2 + ((gy - 0.5) * 1.7) ** 2))
    base = np.array([9, 14, 28]) + gx[..., None] * np.array([4, 14, 22]) + glow[..., None] * np.array([14, 58, 82])
    image = Image.fromarray(base.astype(np.uint8)).convert("RGBA")
    draw = ImageDraw.Draw(image)
    for k in range(-12, 25):
        draw.line([(k * 90 * scale, h), (w * 0.5 + (k - 6) * 22 * scale, h * 0.78)], fill=(40, 78, 104, 90), width=scale)
    for k, y in enumerate([0.79, 0.82, 0.86, 0.91, 0.97]):
        draw.line([(0, h * y), (w, h * y)], fill=(40, 78, 104, 50 + 14 * k), width=scale)
    points = xyz @ rotation(24, -20, -12).T
    depth = 18.0
    projected = [((w * 0.725) + p[0] * 104 * scale * depth / (depth - p[2]), (h * 0.47) - p[1] * 104 * scale * depth / (depth - p[2]), p[2]) for p in points]
    layers = [("atom", projected[i][2], i) for i in range(len(elements))]
    layers += [("bond", (projected[i][2] + projected[j][2]) / 2 - 0.05, (i, j)) for i, j in bonds]
    for kind, _, ref in sorted(layers, key=lambda item: item[1]):
        if kind == "bond":
            i, j = ref
            a, b = np.array(projected[i][:2]), np.array(projected[j][:2])
            direction = (b - a) / np.linalg.norm(b - a)
            ends = [tuple(a + direction * STYLE[elements[i]][0] * 54 * scale), tuple(b - direction * STYLE[elements[j]][0] * 54 * scale)]
            for thickness, shade in [(22, (58, 74, 98)), (14, (116, 138, 166)), (5, (196, 212, 232))]:
                draw.line(ends, fill=shade + (255,), width=thickness * scale)
        else:
            x, y, z = projected[ref]
            size, colour = STYLE[elements[ref]]
            radius = int(size * 68 * scale * depth / (depth - z))
            image.alpha_composite(shaded_sphere(radius, colour), (int(x) - radius, int(y) - radius))
    fonts = Path(matplotlib.get_data_path()) / "fonts" / "ttf"
    face = lambda name, size: ImageFont.truetype(str(fonts / name), size * scale)
    left = 96 * scale
    draw.rectangle([left, 150 * scale, left + 90 * scale, 157 * scale], fill=(72, 202, 228))
    draw.text((left, 184 * scale), "Paracetamol and", font=face("DejaVuSans-Bold.ttf", 70), fill=(244, 248, 252))
    draw.text((left, 272 * scale), "Perioperative Shivering", font=face("DejaVuSans-Bold.ttf", 70), fill=(244, 248, 252))
    draw.text((left, 386 * scale), "Cesarean delivery under spinal anesthesia", font=face("DejaVuSans.ttf", 32), fill=(166, 190, 214))
    draw.text((left, 436 * scale), "Reproducible analysis and data integrity audit", font=face("DejaVuSans.ttf", 32), fill=(166, 190, 214))
    draw.text((left, 510 * scale), "PRELIMINARY AND UNVERIFIED", font=face("DejaVuSans-Bold.ttf", 24), fill=(255, 190, 92))
    draw.text((left, 84 * scale), "R I S E T K U", font=face("DejaVuSans-Bold.ttf", 30), fill=(72, 202, 228))
    image.convert("RGB").resize((width, height), Image.LANCZOS).save(path, optimize=True)


if __name__ == "__main__":
    molecule = paracetamol()
    write_stl(ASSETS / "paracetamol.stl", *molecule)
    banner(ASSETS / "banner.png", *molecule)
    print("atoms", len(molecule[0]), "bonds", len(molecule[2]))
