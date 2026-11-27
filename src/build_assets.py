import random
import struct
from pathlib import Path
import numpy as np

ASSETS = Path(__file__).resolve().parent.parent / "assets"
ASSETS.mkdir(exist_ok=True)
STYLE = {"C": (0.68, "#5b6b86"), "N": (0.66, "#3b76ee"), "O": (0.64, "#ee5546"), "H": (0.42, "#dfe7f2")}
FRAMES, PERIOD = 120, 14
WIDTH, HEIGHT, CENTRE = 1280, 400, (962, 196)
FONT = "font-family=\"'Segoe UI', Helvetica, Arial, sans-serif\""


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
    return elements, xyz - (xyz.max(0) + xyz.min(0)) / 2, bonds


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


def rotation(ax, ay, az):
    a, b, c = np.radians([ax, ay, az])
    rx = np.array([[1, 0, 0], [0, np.cos(a), -np.sin(a)], [0, np.sin(a), np.cos(a)]])
    ry = np.array([[np.cos(b), 0, np.sin(b)], [0, 1, 0], [-np.sin(b), 0, np.cos(b)]])
    rz = np.array([[np.cos(c), -np.sin(c), 0], [np.sin(c), np.cos(c), 0], [0, 0, 1]])
    return rz @ ry @ rx


def project(xyz, turn, scale=47.0, depth=17.0):
    points = xyz @ (rotation(20, 0, -9) @ rotation(0, turn, 0)).T
    factor = depth / (depth - points[:, 2])
    return np.column_stack([CENTRE[0] + points[:, 0] * scale * factor, CENTRE[1] - points[:, 1] * scale * factor, points[:, 2], factor])


def frame(elements, bonds, points):
    items = [(points[i, 2], "atom", i) for i in range(len(elements))]
    items += [((points[i, 2] + points[j, 2]) / 2 - 0.02, "bond", (i, j)) for i, j in bonds]
    out = []
    for _, kind, ref in sorted(items, key=lambda item: item[0]):
        if kind == "atom":
            x, y, _, f = points[ref]
            out.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{STYLE[elements[ref]][0] * 27 * f:.1f}" fill="url(#{elements[ref]})"/>')
        else:
            i, j = ref
            a, b = points[i, :2], points[j, :2]
            length = np.linalg.norm(b - a)
            if length < 1:
                continue
            direction = (b - a) / length
            start = a + direction * min(STYLE[elements[i]][0] * 20, length / 2)
            end = b - direction * min(STYLE[elements[j]][0] * 20, length / 2)
            line = f'x1="{start[0]:.1f}" y1="{start[1]:.1f}" x2="{end[0]:.1f}" y2="{end[1]:.1f}"'
            out.append(f'<line {line} stroke="#45587a" stroke-width="8"/><line {line} stroke="#b9cdea" stroke-width="2.6"/>')
    return "".join(out)


def molecule_layers(elements, xyz, bonds):
    layers = []
    for k in range(FRAMES):
        start, stop = k / FRAMES, (k + 1) / FRAMES
        if k == 0:
            values, times, initial = "visible;hidden", f"0;{stop:.5f}", "visible"
        elif k == FRAMES - 1:
            values, times, initial = "hidden;visible", f"0;{start:.5f}", "hidden"
        else:
            values, times, initial = "hidden;visible;hidden", f"0;{start:.5f};{stop:.5f}", "hidden"
        animate = f'<animate attributeName="visibility" values="{values}" keyTimes="{times}" calcMode="discrete" dur="{PERIOD}s" repeatCount="indefinite"/>'
        layers.append(f'<g visibility="{initial}">{animate}{frame(elements, bonds, project(xyz, 360 * k / FRAMES - 18))}</g>')
    return "".join(layers)


def sphere_gradient(name, colour):
    return (f'<radialGradient id="{name}" cx=".34" cy=".3" r=".75"><stop offset="0" stop-color="#ffffff"/>'
            f'<stop offset=".28" stop-color="{colour}"/><stop offset="1" stop-color="#080d1c"/></radialGradient>')


def pulse_path(y, step=160):
    shape = [(0, 0), (34, 0), (42, -7), (50, 0), (62, 0), (68, 9), (76, -34), (84, 14), (92, 0), (112, 0), (124, -10), (138, 0), (160, 0)]
    points = [(x + k * step, y + dy) for k in range(-1, WIDTH // step + 2) for x, dy in shape]
    return "M" + " L".join(f"{x} {v}" for x, v in points)


def banner(path, elements, xyz, bonds):
    rng = random.Random(20261127)
    stars = "".join(
        f'<circle cx="{rng.randint(10, WIDTH - 10)}" cy="{rng.randint(10, HEIGHT - 10)}" r="{rng.choice([0.9, 1.3, 1.9])}" fill="#7dd3fc">'
        f'<animate attributeName="opacity" values=".1;.9;.1" dur="{rng.uniform(2.4, 5.2):.1f}s" begin="{rng.uniform(0, 3):.1f}s" repeatCount="indefinite"/></circle>'
        for _ in range(46))
    grid = "".join(f'<path d="M0 {y} H{WIDTH}"/>' for y in range(40, HEIGHT, 40))
    orbits = "".join(
        f'<g transform="rotate({tilt} {CENTRE[0]} {CENTRE[1]})"><ellipse cx="{CENTRE[0]}" cy="{CENTRE[1]}" rx="{rx}" ry="{ry}" fill="none" stroke="{colour}" stroke-opacity=".45" stroke-width="1.1" stroke-dasharray="2 9">'
        f'<animate attributeName="stroke-dashoffset" values="0;{shift}" dur="3s" repeatCount="indefinite"/></ellipse>'
        f'<circle r="4" fill="{colour}"><animateMotion dur="{period}s" repeatCount="indefinite" path="M{CENTRE[0] - rx} {CENTRE[1]} a{rx} {ry} 0 1 0 {2 * rx} 0 a{rx} {ry} 0 1 0 {-2 * rx} 0"/></circle></g>'
        for tilt, rx, ry, colour, period, shift in [(-14, 262, 92, "#22d3ee", 9, 110), (24, 236, 66, "#a78bfa", 12, -110)])
    chips, x = [], 90
    for label, width in [("Open label trial", 138), ("Reproducible notebook", 188), ("Data integrity audit", 170)]:
        chips.append(f'<rect x="{x}" y="298" width="{width}" height="32" rx="16" fill="#0b1b30" fill-opacity=".7" stroke="#22d3ee" stroke-opacity=".55"/>'
                     f'<text x="{x + width / 2}" y="319" text-anchor="middle" font-size="14" fill="#e0f2fe" {FONT}>{label}</text>')
        x += width + 14
    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {WIDTH} {HEIGHT}" width="{WIDTH}" height="{HEIGHT}" role="img" aria-labelledby="t d">
<title id="t">Risetku paracetamol and shivering study</title><desc id="d">Animated three dimensional paracetamol molecule rotating beside the study title.</desc>
<defs>
<linearGradient id="bg" x1="0" y1="0" x2="1" y2="1"><stop offset="0" stop-color="#060a14"/><stop offset=".55" stop-color="#0b1330"/><stop offset="1" stop-color="#130c2e"/></linearGradient>
<radialGradient id="halo" cx=".75" cy=".5" r=".42"><stop offset="0" stop-color="#22d3ee" stop-opacity=".3"/><stop offset="1" stop-color="#22d3ee" stop-opacity="0"/></radialGradient>
<radialGradient id="halo2" cx=".15" cy=".15" r=".55"><stop offset="0" stop-color="#7c3aed" stop-opacity=".26"/><stop offset="1" stop-color="#7c3aed" stop-opacity="0"/></radialGradient>
<linearGradient id="beam" x1="0" x2="1"><stop offset="0" stop-color="#22d3ee" stop-opacity="0"/><stop offset=".5" stop-color="#a5f3fc" stop-opacity=".95"/><stop offset="1" stop-color="#22d3ee" stop-opacity="0"/></linearGradient>
<linearGradient id="title" x1="0" x2="1"><stop offset="0" stop-color="#ffffff"/><stop offset=".5" stop-color="#a5f3fc"/><stop offset="1" stop-color="#c4b5fd"/>
<animate attributeName="x1" values="-.6;.6;-.6" dur="7s" repeatCount="indefinite"/><animate attributeName="x2" values=".4;1.6;.4" dur="7s" repeatCount="indefinite"/></linearGradient>
{"".join(sphere_gradient(k, v[1]) for k, v in STYLE.items())}
<filter id="soft" x="-50%" y="-200%" width="200%" height="500%"><feGaussianBlur stdDeviation="9"/></filter>
</defs>
<rect width="{WIDTH}" height="{HEIGHT}" fill="url(#bg)"/><rect width="{WIDTH}" height="{HEIGHT}" fill="url(#halo)"/><rect width="{WIDTH}" height="{HEIGHT}" fill="url(#halo2)"/>
<g opacity=".13" stroke="#67e8f9" stroke-width=".6">{grid}<animateTransform attributeName="transform" type="translate" values="0 0;0 40" dur="6s" repeatCount="indefinite"/></g>
{stars}
<path d="{pulse_path(366)}" fill="none" stroke="#22d3ee" stroke-opacity=".5" stroke-width="1.6"><animateTransform attributeName="transform" type="translate" values="0 0;-160 0" dur="2.6s" repeatCount="indefinite"/></path>
<ellipse cx="{CENTRE[0]}" cy="338" rx="196" ry="16" fill="#22d3ee" opacity=".28" filter="url(#soft)"><animate attributeName="rx" values="196;160;196" dur="5s" repeatCount="indefinite"/></ellipse>
{orbits}
<g><animateTransform attributeName="transform" type="translate" values="0 -7;0 7;0 -7" dur="5s" repeatCount="indefinite"/>{molecule_layers(elements, xyz, bonds)}</g>
<rect y="0" x="720" width="480" height="2" fill="url(#beam)" opacity=".7"><animate attributeName="y" values="56;336;56" dur="6.5s" repeatCount="indefinite"/></rect>
<text x="92" y="84" font-size="15" letter-spacing="7" fill="#67e8f9" {FONT}>RISETKU CLINICAL RESEARCH</text>
<text x="88" y="154" font-size="60" font-weight="800" letter-spacing="-1" fill="url(#title)" {FONT}>Paracetamol and</text>
<text x="88" y="218" font-size="60" font-weight="800" letter-spacing="-1" fill="url(#title)" {FONT}>Shivering</text>
<text x="92" y="254" font-size="19" fill="#dbeafe" {FONT}>Cesarean delivery under spinal anesthesia</text>
<text x="92" y="279" font-size="19" fill="#dbeafe" {FONT}>A reproducible analysis with data integrity screening</text>
{"".join(chips)}
</svg>'''
    Path(path).write_text(svg, encoding="utf-8")


if __name__ == "__main__":
    molecule = paracetamol()
    write_stl(ASSETS / "paracetamol.stl", *molecule)
    banner(ASSETS / "banner.svg", *molecule)
    print("atoms", len(molecule[0]), "bonds", len(molecule[2]), "banner kB", (ASSETS / "banner.svg").stat().st_size // 1024)
