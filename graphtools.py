from math import sqrt, sin, cos, pi

class CubicVertex:
	def __init__(self, graph, branch):
		self.graph = graph
		self.name = branch[0]
		self.neighbours_names = list(branch[1:])
		self.obey = 1

	def switch_obedience(self):
		self.obey = -self.obey

	def set_to_disobey(self):
		self.obey = -1

	def set_to_obey(self):
		self.obey = 1

	def __repr__(self):
		r = f"Vertex '{self.name}' of CubicGraph #{self.graph.uid} (neighbours={'|'.join([self[i] for i in [0, 1, 2]])})"
		return r

	def __getitem__(self, n):
		# print(f"Looking for neighbour {repr(n)} of vertex '{self.name}'.")
		try:
			return self.neighbours_names.index(n)
		except ValueError:
			return self.neighbours_names[n%3]

	def next(self, n, spin=1):
		return self[(self[n]+spin)%3]
			
class CubicGraph:

	uid = 0
	
	def __init__(self, *direct_branches):
		self.uid = type(self).uid
		type(self).uid += 1
		self.vertices = {b[0]:CubicVertex(self, b) for b in direct_branches}
		self.ordered_edges = [[u, v] for u in self.vertices for v in self[u].neighbours_names]
		self.edges = [[u, v] for u, v in self.ordered_edges if u < v]

	def __len__(self):
		return len(self.vertices)

	def __repr__(self):
		r = f"CyclicGraph #{self.uid}, VEV*={(len(self), 3*len(self)//2, len(self)//2+2)}, vertices = {'|'.join(self.vertices.keys())}"
		return r

	def __getitem__(self, v):
		# (f"Looking for {v} in graph.")
		return self.vertices[v]

	def set_to_obey(self, *vertices):
		for v in vertices:
			self[v].set_to_obey()
			
	def set_to_disobey(self, *vertices):
		for v in vertices:
			self[v].set_to_disobey()
			
	def set_obedience(self, *disobeying):
		for v in self.vertices:
			if v in disobeying:
				print(f"'{v}' disobeys")
				self[v].set_to_disobey()
			else:
				print(f"'{v}' obeys")
				self[v].set_to_obey()

	def find_right_cycles(self):
		cycles = []
		for v in self.vertices:
			for u in self[v].neighbours_names:
				if all(v+u not in c+c[0] for c in cycles):
					w = SignedWalk(self, v, u)
					cycles.append(w.vertices)
		return cycles

	def dual(self):
		dual_vertices = {i:c for i, c in enumerate(self.find_right_cycles())}
		edge_map = {v+u:[] for v in self.vertices for u in self[v].neighbours_names if v<u}
		for edge in edge_map:
			for i, c in dual_vertices.items():
				if edge in c+c[0] or edge[::-1] in c+c[0]:
					edge_map[edge].append(i)
		return dual_vertices, edge_map
			
class Vector:

	epsilon = 0.0000001

	@classmethod
	def copy(cls, v):
		return cls(v.x, v.y)

	@classmethod
	def polar(cls, r, t):
		return cls(r*cos(t/180*pi), r*sin(t/180*pi))
	
	def __init__(self, x, y):
		self.x = x
		self.y = y

	def __repr__(self):
		return f"Vector({self.x, self.y})"

	def __abs__(self):
		return sqrt(self.x*self.x+self.y*self.y)

	def __add__(self, other):
		return Vector(self.x+other.x, self.y+other.y)

	def __sub__(self, other):
		return Vector(self.x-other.x, self.y-other.y)

	def __mul__(self, other):
		return self.x*other.x + self.y*other.y

	def __rmul__(self, r):
		return Vector(self.x*r, self.y*r)

	def __truediv__(self, r):
		return Vector(self.x/r, self.y/r)

	def __floordiv__(self, other):
		d = abs(self.x*other.y - self.y*other.x)/(abs(self)*abs(other))
		return d<self.epsilon

	def __imul__(self, r):
		self.x *= r
		self.y *= r

	def __itruediv__(self, r):
		self.x /= r
		self.y /= r
		return self

	def normalize(self):
		self /= abs(self)
		return self

	def normal(self):
		n = Vector(-self.y, self.x)
		n.normalize()
		return n

	def print_coords(self, precision=4):
		return f"({float(self.x):.4}, {float(self.y):.4})"

class TikZGraph:

	
	def __init__(self, graph, *, cartesian={}, polar={}, vector={}, scale=1, shift=.05, edge_delta = .2, center_delta = .1):
		self.graph = graph
		self.shift = shift
		self.edge_delta = edge_delta
		self.center_delta = center_delta
		c_dict = {v:scale*Vector(*c) for v, c in cartesian.items()}
		p_dict = {v:scale*Vector.polar(*c) for v, c in polar.items()}
		v_dict = {v:scale*Vector.copy(c) for v, c in vector.items()}
		self.vertex_vectors = {**c_dict, **p_dict, **v_dict}
		self.vertex_anchors = {}
		self.oriented_edge_anchors = {}

		# generating anchor points : for every vertex and every edge
		#
		# /----------------7-----------------\
		#   0   1        2 3 4          5   6
		for nb in self.graph.vertices:
			self.oriented_edge_anchors[nb] = {}
			for i in [0, 1, 2]:
				na = self.graph[nb][i-1]
				nc = self.graph[nb][i]
				nd = self.graph.vertices[nc].next(nb)
				a, b, c, d = (self.vertex_vectors[w] for w in (na, nb, nc, nd))

				sa, *coefs_ab = self.shifted_equation(a, b)
				sb, *coefs_bc = self.shifted_equation(b, c)
				sc, *coefs_cd = self.shifted_equation(c, d)

				if (b-a)//(c-b):
					anchor0 = sb
				else:
					anchor0 = self.intersection_of(coefs_ab, coefs_bc)

				if (c-b)//(d-c):
					anchor6 = sc
				else:
					anchor6 = self.intersection_of(coefs_bc, coefs_cd)

				director = (anchor6 - anchor0).normalize()
				anchor3 = (anchor0 + anchor6)/2

				anchor1 = anchor0 + self.edge_delta*director
				anchor2 = anchor3 - self.center_delta*director
				anchor4 = anchor3 + self.center_delta*director
				anchor5 = anchor6 - self.edge_delta*director

				self.oriented_edge_anchors[nb][nc] = [anchor0, anchor1, anchor2, anchor3, anchor4, anchor5, anchor6, (b+c)/2]

	def generate_tikz_graph(self, *, insert=[], options="", environment=False, centered=True):
		s = ["% ------ Graph -----------"]
		if environment:
			if centered:
				s.append("\\begin{center}")
			s.append(f"\\begin{{tikzpicture}}[{options}]")
		for v in self.graph.vertices:
			s.append(f"\coordinate (node{v}) at {self.vertex_vectors[v].print_coords()};")
		s.append("\\draw ")
		for a, b in self.graph.edges:
			s.append(f" (node{a}) -- (node{b})")
		s.append(";")
		s.extend(insert)
		if environment:
			s.append("\end{tikzpicture}")
			if centered:
				s.append("\end{center}")
			
		return "\n".join(s)

	def generate_tikz_vertex_labels(self):
		s = ["% ------ Vertex labels ---"]
		for n, v in self.graph.vertices.items():
			s.append(f"\\node[vertex, {'dis' if v.obey==-1 else ''}obey] at (node{n}) {{{n}}};")
		return "\n".join(s)

	def generate_tikz_walk(self, walk, *, arrows=True, options=""):
		a, b, spin = walk.steps[-1]
		if spin == 1:
			p = self.oriented_edge_anchors[a][b][5]
		else:
			p = self.oriented_edge_anchors[b][a][1]
		s = ["% ------ Walk ------------", f"\\draw[{options}] {p.print_coords()} "]
		for c, d, next_spin in walk.steps:
			da = self.oriented_edge_anchors[c][d]
			db = self.oriented_edge_anchors[d][c]
			if spin == 1:
				if spin == next_spin:
					code = f".. controls {da[0].print_coords()} .. {da[1].print_coords()} -- "
					if arrows:
						code += f"pic[pos=0]{{arrow}} "
					code += f"{da[5].print_coords()}"
					s.append(code)
				else:
					code = f".. controls {da[0].print_coords()} .. {da[1].print_coords()} -- "
					if arrows:
						code += f"pic[pos=0]{{arrow}} "
					code += f"{da[2].print_coords()}" \
						f" .. controls {da[3].print_coords()} and {db[3].print_coords()} .. {db[2].print_coords()}" \
						f" -- {db[1].print_coords()}"
					s.append(code)
			else:
				if spin == next_spin:
					code = f".. controls {db[6].print_coords()} .. {db[5].print_coords()} -- "
					if arrows:
						code += f"pic[pos=0]{{arrow}} "
					code += f"{db[1].print_coords()}"
					s.append(code)
				else:
					code = f".. controls {db[6].print_coords()} .. {db[5].print_coords()} -- "
					if arrows:
						code += f"pic[pos=0]{{arrow}} "
					code += f"{db[4].print_coords()}" \
						f" .. controls {db[3].print_coords()} and {da[3].print_coords()} .. {da[4].print_coords()}" \
						f" -- {da[5].print_coords()}"
					s.append(code)
			a, b, spin = c, d, next_spin
		s.append(";")
		return "\n".join(s)

	def intersection_of(self, line_a, line_b):
		ax, ay, ad = line_a
		bx, by, bd = line_b
		
		#  ax * x + ay * y + ad = 0
		#  bx * x + by * y + bd = 0
		#
		# (by*ax - ay*bx) x + (by*ay - ay*by) y + by*ad - ay*bd = 0 
		
		x = (ay*bd - by*ad)/(by*ax - ay*bx)
		y = (ax*bd - bx*ad)/(bx*ay - ax*by)

		return Vector(x, y)
				
	def shifted_equation(self,a, b):
		
		director = b-a
		normal = director.normal()
		director.normalize()

		sa = a - self.shift*normal
		sb = b - self.shift*normal

		# line coefs
		cx = (sa.y - sb.y)
		cy = (sb.x - sa.x)
		d = -cx*sa.x -cy*sa.y

		return sa, cx, cy, d
		
	def __getitem__(self, request):
		try:
			return self.vertex_vectors[request]
		except KeyError:
			print("Beep!", request)

	def compute_dual(self, dual_vertices, dual_edges, outer_cell = None):
		dv_data = {}
		for i, vertices in dual_vertices.items():
			if i!=outer_cell:
				x = sum(self.vertex_vectors[c].x for c in vertices)/len(vertices)
				y = sum(self.vertex_vectors[c].y for c in vertices)/len(vertices)
				dv_data[i] = Vector(x, y)
				print(f"\\coordinate (dual{i}) at {(x,y)};")

		outer_edges = []

		if outer_cell is not None:
			vertices = dual_vertices[outer_cell]
			for k in range(len(vertices)):
				s, e = vertices[k-1], vertices[k]
				if s > e:
					s, e = e, s
				a, b = dual_edges[s+e]
				if a == outer_cell:
					a = b
				outer_edges.append(a)
 
				se = self[e]-self[s]
				m = dv_data[a]
				sm = m - self[e]
				smp = (2*(sm*se)*se)/(se*se) - sm
				mp = smp + self[s]

				print(f"\\coordinate (dual out{a}) at {mp.print_coords()};")
					 
		for edge in dual_edges:
			if not outer_cell in dual_edges[edge]:
				a, b = dual_edges[edge]
				print(f"\\draw[dual edge] (dual{a}) -- (dual{b});")

		for edge in outer_edges:
			print(f"\\draw[dual edge] (dual{edge}) -- (dual out{edge});")


	vertex_type_labels = {
		"A1" : "A\\textsubscript{1}",
		"A2" : "A\\textsubscript{2}",
		"B1" : "B\\textsubscript{1}",
		"B2" : "B\\textsubscript{2}",
		"D1" : "D\\textsubscript{1}",
		"D2" : "D\\textsubscript{2}",
		"d" : "$\overline{\mbox{D}}$",
		"*" : ""
	}

	def generate_vertex_types(self, walk, options=""):
		s = ["% ------ Vertex types ---"]
		types = walk.vertex_types
		for n, v in self.graph.vertices.items():
			s.append(f"\\node[vertex type, {options}] at (node{n}) {{{self.vertex_type_labels[types[n]]}}};")
		return "\n".join(s)
		
		

class SignedWalk:

	@classmethod
	def alternate(self, graph, start, aim, spin=1):
		return SignedWalk(graph, start, aim, [spin, -spin])

	def __init__(self, graph, start, aim, signature = [1]):
		self.graph = graph
		self.steps = []
		self.raw_steps = []
		spin = signature[0]
		count = 0
		l = len(signature)
		while not (start, aim, spin) in self.raw_steps:
			self.raw_steps.append((start, aim, spin))
			spin *= graph[aim].obey
			self.steps.append((start, aim, spin))
			
			count += 1
			start, aim, spin = aim, graph[aim].next(start, spin), signature[count%l]

	@property
	def vertices(self):
		try:
			return self.__vertices
		except AttributeError:
			self.__vertices = "".join(s[0] for s in self.steps)
			return self.__vertices

	def summary(self):
		nloops = 0
		state = {v:0 for v in self.graph.vertices}
		freq = {v:{0:{u:0 for u in self.graph.vertices if u!=v}} for v in self.graph.vertices}
		for v, n, s in self.steps:
			state[v] += 1
			freq[v][state[v]] = {u:0 for u in self.graph.vertices if u != v}
			for u in self.graph.vertices:
				if u != v:
					freq[u][state[u]][v] += 1
				
		for v in self.graph.vertices:
			if 3 in freq[v]:
				for u in freq[v][3]:
					freq[v][0][u] += freq[v][3][u]
				del freq[v][3]
		
		return freq

	@property
	def vertex_types(self):
		try:
			return self.__vertex_types
		except AttributeError:
			self.__vertex_types = self.__compute_vertex_types()
			return self.__vertex_types


	def __compute_vertex_types(self):
		type_of = {}
		vs = self.vertices
		n = len(vs)
		for v in self.graph.vertices:
			encounters = []
			i = vs.find(v)
			while i >-1 :
				encounters.append(vs[i-1] + vs[i] + vs[(i+1)%n])
				i = vs.find(v, i+1)
			if len(encounters) == 0:
				type_of[v] = "*"
			elif len(encounters) == 1:
				type_of[v] = "d"
			elif len(encounters) == 2:
				if encounters[0][2] == encounters[1][0]:
					type_of[v] = "D1"
				else:
					type_of[v] = "D2"
			else:
				if all(encounters[i-1][2] == encounters[i][0] for i in range(3)):
					type_of[v] = "A1"
				elif all(encounters[i-1][2] == encounters[(i+1)%3][0] for i in range(3)):
					type_of[v] = "A2"
				elif any(encounters[i-1][2] == encounters[i][0] for i in range(3)):
					type_of[v] = "B2"
				else:
					type_of[v] = "B1"
		
		return type_of

