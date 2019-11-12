#include <iostream>
#include <algorithm>
#include <string.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <cassert>
#include  <cmath>
#include <GL/glut.h>
using namespace std;

// -------------- parameters -------------- 
const int width = 800;
const int height = 800;
float timer = 0;

int mousehitX, mousehitY;
int mouseX, mouseY;
int button_kind = -1;
bool mousePressed = false;

bool hasDrawn = false;

int frameBufferR[width][height] = { 0 };
int frameBufferG[width][height] = { 0 };
int frameBufferB[width][height] = { 0 };


// -------------- Mouse IO -------------------

enum BtnState {
	IDLE,
	HOVER,
	PRESSED
};
class Button {
public:
	BtnState state;
	int p1x, p1y, p2x, p2y;
	int c1r, c1g, c1b;
	int c2r, c2g, c2b;
	int c3r, c3g, c3b;
	const char* btn_str;

	Button(int _p1x, int _p1y, int _p2x, int _p2y, const char* _btn_str) {
		p1x = _p1x;
		p1y = _p1y;
		p2x = _p2x;
		p2y = _p2y;

		btn_str = _btn_str;

		c1r = 96;
		c1b = 120;
		c1g = 130;
		c2r = 40;
		c2g = 150;
		c2b = 250;
		c3r = 150;
		c3g = 230;
		c3b = 255;
	}

	void drawText(const char* str, int x, int y) {
		int n = strlen(str);
		glRasterPos2i(x, y);
		for (int i = 0; i < n; i++)
			glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, *(str + i));

	}

	void updateBtn() {
		if (mouseX > p1x&& mouseX < p2x && mouseY > p1y&& mouseY < p2y) {
			if (mousePressed)
				state = PRESSED;
			else
				state = HOVER;
		}
		else {
			state = IDLE;
		}
	}
	void drawBtn() {
		updateBtn();

		switch (state) {
		case IDLE:
			glColor3f(c1r / 255., c1g / 255., c1b / 255.);
			break;
		case HOVER:
			glColor3f(c2r / 255., c2g / 255., c2b / 255.);
			break;
		case PRESSED:
			glColor3f(c3r / 255., c3g / 255., c3b / 255.);
			break;
		default:
			glColor3f(c1r / 255., c1g / 255., c1b / 255.);
			break;
		}
		glRecti(p1x, p1y, p2x, p2y);

		drawText(btn_str, p1x, p1y - 15);
	}
};
void mouse_hit(int button, int state, int x, int y) {
	button_kind = button;

	switch (button) {
	case GLUT_LEFT_BUTTON:
		if (state == GLUT_DOWN) {
			mousehitX = x;
			mousehitY = height - y;
			mousePressed = true;

		}
		else if (state == GLUT_UP) {
			mousePressed = false;

			cout << ">>>>>>>>>>>>  Processing Data...  <<<<<<<<<<<" << endl;

			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
			glClearColor(0, 0, 0, 1);
			memset(frameBufferR, 0, sizeof(frameBufferR));
			memset(frameBufferG, 0, sizeof(frameBufferG));
			memset(frameBufferB, 0, sizeof(frameBufferB));
			hasDrawn = false;
		}
		break;

	case GLUT_RIGHT_BUTTON:
		if (state == GLUT_DOWN) {
			mousehitX = x;
			mousehitY = height - y;
		}
		break;

	default:
		break;
	}

}
void mouse_move(int x, int y) {
	mouseX = x;
	mouseY = height - y;
}
vector<Button*> shader_btns;
vector<Button*> model_btns;
int shaderActiveBtn = 0;
int modelActiveBtn = 2;

void initUI() {
	//  --------- Add Elements in UI ---------------
	Button* phong = new Button(width - 160, 50, width - 80, 75, "Phong Shading");
	Button* gourand = new Button(width - 160, 100, width - 80, 125, "Gourand Shading");
	Button* halftone = new Button(width - 160, 150, width - 80, 175, "Half-Tone Shading");
	Button* toon = new Button(width - 160, 200, width - 80, 225, "Toon Shading");
	shader_btns.push_back(phong);
	shader_btns.push_back(gourand);
	shader_btns.push_back(halftone);
	shader_btns.push_back(toon);

	Button* bunny = new Button(width - 280, 50, width - 200, 75, "Bunny");
	Button* cube = new Button(width - 280, 100, width - 200, 125, "Cube Icosahedron");
	Button* ven = new Button(width - 280, 150, width - 200, 175, "Beethoven");
	model_btns.push_back(bunny);
	model_btns.push_back(cube);
	model_btns.push_back(ven);
}


// -------------- Basic Structs -------------- 
struct Color {
	Color() : r(0), g(0), b(0) {}
	Color(int gray) : r(gray), g(gray), b(gray) {}
	Color(int _r, int _g, int _b) : r(_r), g(_g), b(_b) {}
	int r, g, b;
};
struct Colorf {
	Colorf() : r(0), g(0), b(0) {}
	Colorf(float gray) : r(gray), g(gray), b(gray) {}
	Colorf(float _r, float _g, float _b) : r(_r), g(_g), b(_b) {}

	Colorf operator *(const float f) { return Colorf(r * f, g * f, b * f); }
	Colorf operator /(const float f) { return Colorf(r / f, g / f, b / f); }

	float r, g, b;
};

struct Vertex {
	Vertex() : x(0), y(0), z(0) {}
	Vertex(float _x, float _y) : x(_x), y(_y), z(0) {}
	Vertex(float _x, float _y, float _z) : x(_x), y(_y), z(_z) {}

	Vertex operator +(const Vertex& v) const { return Vertex(x + v.x, y + v.y, z + v.z); }
	Vertex operator -(const Vertex& v) const { return Vertex(x - v.x, y - v.y, z - v.z); }
	Vertex operator *(float f) const { return Vertex(x * f, y * f, z * f); }
	Vertex operator /(float f) const { return Vertex(x / f, y / f, z / f); }
	float operator [](const int i) { assert(i < 3); return i <= 0 ? x : (1 == i ? y : z); }
	//const float operator[](const int i) const { assert(i < 3); return i <= 0 ? x : (1 == i ? y : z); }

	float norm() const { return sqrt(x * x + y * y + z * z); }
	//Dot Product
	float operator *(const Vertex& v) { return x * v.x + y * v.y + z * v.z; }
	//Cross Product
	Vertex operator ^(const Vertex& v) const { return Vertex(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x); }
	// Nornalize
	Vertex& normalize(float l = 1) { *this = (*this) * (l / norm()); return *this; }

	float x, y, z;
};

struct Face {
	Face() : f1(0), f2(0), f3(0) {}
	Face(int _f1, int _f2, int _f3) : f1(_f1), f2(_f2), f3(_f3) {}
	int operator [](const int i) { assert(i < 3);  return i <= 0 ? f1 : (1 == i ? f2 : f3); }

	int f1, f2, f3;
};
struct Geometry {
	Geometry() : verts(vector<Vertex>()), norms(vector<Vertex>()), colors(vector<Color>()), faces(vector<Face>()), specs(vector<int>()) {}
	Geometry(vector<Vertex> _verts, vector<Color> _colors, vector<Face> _faces, vector<int> _specs) : verts(_verts), norms(vector<Vertex>()), colors(_colors), faces(_faces), specs(_specs) {}
	vector<Vertex> verts;
	vector<Vertex> norms;
	vector<Color> colors;
	vector<Face> faces;
	vector<int> specs;

	// Calculate the Average Z depth of one face
	float getXYZAvg(Face face, int axis) {
		float axisSum = 0;
		for (int j = 0; j < 3; j++) {
			Vertex world_coords = verts[face[j] - 1];
			if (axis == 0)
				axisSum += world_coords.z;
			if (axis == 1)
				axisSum += world_coords.y;
			if (axis == 2)
				axisSum += world_coords.x;
		}
		return axisSum / 3.0;
	}
	// Sort Faces by depth in order to draw like a "Painter"~ lol
	void sortFaces(int axis) {
		cout << "Sort Good!" << endl;
		for (int i = 0; i < faces.size() - 1; i++) {
			for (int j = 0; j < faces.size() - 1 - i; j++) {
				if (getXYZAvg(faces[j], axis) > getXYZAvg(faces[j + 1], axis)) {
					swap(faces[j], faces[j + 1]);
				}
			}
		}
	}
	// calculate normals of each vertex by averaging all face normals connecting to it.
	void calcNormals() {
		cout << "Calc Norm Good!" << endl;
		for (int i = 0; i < verts.size(); i++) {
			vector<int> inWhichFaces;
			for (int j = 0; j < faces.size(); j++) {
				for (int k = 0; k < 3; k++) {
					if (faces[j][k] - 1 == i) {
						inWhichFaces.push_back(j);
					}
				}
			}
			vector<Vertex> facesNormals;
			for (int m = 0; m < inWhichFaces.size(); m++) {
				Face thisFace = faces[inWhichFaces[m]];
				Vertex v1 = verts[thisFace[0] - 1];
				Vertex v2 = verts[thisFace[1] - 1];
				Vertex v3 = verts[thisFace[2] - 1];
				Vertex N;
				N = ((v2 - v1) ^ (v3 - v1));
				N.normalize();
				facesNormals.push_back(N);
			}
			Vertex vAvg;
			for (int m = 0; m < facesNormals.size(); m++) {
				vAvg = vAvg + facesNormals[m];
			}
			vAvg / facesNormals.size();
			norms.push_back(vAvg);
		}
	}
	// Get the Index of normal by face index
	int getNormIndex(int iface, int nvert) {
		return faces[iface][nvert];
	}
	// Get the normal vector by face index
	Vertex getNormal(int iface, int nvert) {
		return norms[faces[iface][nvert] - 1];
	}
};

// --------------- Model --------------------
class Model {
public:
	// Real Container
	vector<Geometry> geos;

	~Model();
	Model() {

	}
	Model(const char* _filename) {
		loadFile(_filename);
	}
	Model(const Model& m) {
		geos = m.geos;
	}

	vector<string>& string_split(const string& str, char delim, vector<string>& elems, bool skip_empty = true) {
		istringstream iss(str);
		for (string item; getline(iss, item, delim); )
			if (skip_empty && item.empty()) continue;
			else elems.push_back(item);
		return elems;
	}
	void loadFile(const char* filename) {
		// Temp Container
		vector<Vertex> verts;
		vector<Vertex> norms;
		vector<Color> colors;
		vector<Face> faces;
		vector<int> specs;


		ifstream in(filename);
		//in.open(filename, ifstream::in);
		if (in.fail())
			return;
		string line;

		getline(in, line); // Get Geometry Nums
		int geoNums = stoi(line);
		int geoIndex = geoNums;
		cout << geoNums << endl;

		getline(in, line); // skip a blank line

		getline(in, line);
		int vertsNum = stoi(line);
		int colorNum = vertsNum;
		cout << vertsNum << endl;
		int faceNum = -1;
		int specsNum = -1;
		while (!in.eof() && geoIndex > 0) {
			if (vertsNum > 0) {
				//cout << "read vert" << vertsNum << endl;
				getline(in, line);
				vector<string> result;
				string_split(line, ' ', result);
				verts.push_back(Vertex(stof(result[0]), stof(result[1]), stof(result[2])));
				vertsNum--;
			}
			else if (colorNum > 0) {
				//cout << "read color" << colorNum << endl;
				getline(in, line);
				vector<string> result;
				string_split(line, ' ', result);
				colors.push_back(Color(stoi(result[0]), stoi(result[1]), stoi(result[2])));
				colorNum--;
			}
			else if (colorNum == 0 && faceNum == -1 && specsNum == -1) {
				getline(in, line);
				faceNum = stoi(line);
				specsNum = stoi(line);
			}
			else if (faceNum > 0) {
				//cout << "Read Face: " << faceNum << endl;
				getline(in, line);
				vector<string> result;
				string_split(line, ' ', result);
				faces.push_back(Face(stoi(result[0]), stoi(result[1]), stoi(result[2])));
				faceNum--;
			}
			else if (specsNum > 0) {
				//cout << "Read Spec: " << specsNum << endl;
				getline(in, line);
				vector<string> result;
				string_split(line, ' ', result);
				specs.push_back(stoi(result[0]));
				specsNum--;
			}
			else {
				geoIndex--;
				cout << "Face Num: " << faces.size() << " Specs Num: " << specs.size() << endl;
				geos.push_back(Geometry(verts, colors, faces, specs));
				if (geoIndex > 0) {
					getline(in, line);	// skip a blank line
					getline(in, line);
					vertsNum = stoi(line);
					colorNum = vertsNum;
					faceNum = -1;
					specsNum = -1;
					verts.clear();
					colors.clear();
					faces.clear();
					specs.clear();
				}
			}
		}
	}

	vector<Geometry> getGeo() { return geos; };
};



// -------------- Verts and Faces (Deprecated) -------------------
vector<Vertex> verts;
vector<Vertex> norms;
vector<Color> colors;
vector<Face> faces;
vector<int> specs;







// ------------- Tool Functions (Deprecated) ----------------
float lerp(float _a, float _b, float _t) {
	return _a + (_b - _a) * _t;
}
vector<string>& string_split(const string& str, char delim, vector<string>& elems, bool skip_empty = true) {
	istringstream iss(str);
	for (string item; getline(iss, item, delim); )
		if (skip_empty && item.empty()) continue;
		else elems.push_back(item);
	return elems;
}
// Standalone cross product function  z = x ^ y
Vertex crossProduct(const Vertex x, const Vertex y) {
	Vertex z;
	float m1, m2, m3;
	m1 = x.y * y.z - x.z * y.y;
	m2 = x.z * y.x - x.x * y.z;
	m3 = x.x * y.y - x.y * y.x;
	z.x = m1;
	z.y = m2;
	z.z = m3;
	return z;
}

// calculate normals of each vertex by averaging all face normals connecting to it.
void calcNormals() {
	for (int i = 0; i < verts.size(); i++) {
		vector<int> inWhichFaces;
		for (int j = 0; j < faces.size(); j++) {
			for (int k = 0; k < 3; k++) {
				if (faces[j][k] - 1 == i) {
					inWhichFaces.push_back(j);
				}
			}
		}
		vector<Vertex> facesNormals;
		for (int m = 0; m < inWhichFaces.size(); m++) {
			Face thisFace = faces[inWhichFaces[m]];
			Vertex v1 = verts[thisFace[0] - 1];
			Vertex v2 = verts[thisFace[1] - 1];
			Vertex v3 = verts[thisFace[2] - 1];
			Vertex N;
			N = crossProduct((v2 - v1), (v3 - v1));
			N.normalize();
			facesNormals.push_back(N);
		}
		Vertex vAvg;
		for (int m = 0; m < facesNormals.size(); m++) {
			vAvg = vAvg + facesNormals[m];
		}
		vAvg / facesNormals.size();
		norms.push_back(vAvg);
	}
}
// Get the Index of normal by face index
int getNormIndex(int iface, int nvert) {
	return faces[iface][nvert];
}
// Get the normal vector by face index
Vertex getNormal(int iface, int nvert) {
	return norms[faces[iface][nvert] - 1];
}
// Calculate the Average Z depth of one face
float getZAvg(Face face) {
	float zSum = 0;
	for (int j = 0; j < 3; j++) {
		Vertex world_coords = verts[face[j] - 1];
		zSum += world_coords.z;
	}
	return zSum / 3.0;
}
// Sort Faces by depth in order to draw like a "Painter"~ lol
void sortFaces() {
	for (int i = 0; i < faces.size() - 1; i++) {
		for (int j = 0; j < faces.size() - 1 - i; j++) {
			if (getZAvg(faces[j]) > getZAvg(faces[j + 1])) {
				swap(faces[j], faces[j + 1]);
			}
		}
	}
}


//Vertex world2screen(Vertex v) {
//	return Vertex(int((v.x - 0.5) * width * 2 + width / 2), int((v.y - 0.5) * height * 2 + height / 2), int(v.z));
//}
Vertex world2screen(Vertex v) {
	return Vertex(int(v.x * width), int(v.y * height), int(v.z));
}

// --------------- Shader Tool Functions --------------
float step(float threshould, float val) {
	return val >= threshould;
}
float frac(float v) {
	return v - floor(v);
}
float saturate(float x) {
	return max(0.0f, min(1.0f, x));
}
float smoothstep(float a, float b, float x) {
	float t = saturate((x - a) / (b - a));
	return t * t * (3.0f - (2.0f * t));
}
float circle(Vertex _st, float _radius) {
	Vertex pos = Vertex(0.5f, 0.5f) - _st;
	return smoothstep(1.0 - _radius, 1.0 - _radius + _radius * 0.2, 1. - (pos * pos) * 3.14);
}
// --------------- Shaders -----------------------
struct IShader {
	virtual ~IShader();
	virtual Vertex vertexShader(int iface, int nvert) = 0;
	virtual Colorf fragmentShader(Vertex uv, Vertex barycentric) = 0;
	virtual void setG(Geometry _g) = 0;
};
struct GourandShader : public IShader {
	float varying_intensity[3];
	int varying_vert;
	Vertex light_dir;
	Geometry g;

	virtual Vertex vertexShader(int iface, int nvert) {
		varying_vert = nvert;
		Vertex world_coord = g.verts[g.faces[iface][nvert] - 1];
		varying_intensity[nvert] = g.getNormal(iface, nvert) * light_dir;
		return world_coord;
	}
	virtual Colorf fragmentShader(Vertex uv, Vertex barycentric) {
		Vertex vi = Vertex(varying_intensity[0], varying_intensity[1], varying_intensity[2]);
		float intensity = vi * barycentric;
		intensity *= 0.13;

		Colorf baseColor = Colorf(g.colors[varying_vert].r / 255.0, g.colors[varying_vert].g / 255.0, g.colors[varying_vert].b / 255.0);
		Colorf result = baseColor * intensity;
		return result;
	}
	virtual void setG(Geometry _g) {
		g = _g;
	}
};
struct ToonShader : public IShader {
	float varying_intensity[3];
	int varying_vert;
	Vertex light_dir;
	Geometry g;

	virtual Vertex vertexShader(int iface, int nvert) {
		varying_vert = nvert;
		Vertex world_coord = g.verts[g.faces[iface][nvert] - 1];
		varying_intensity[nvert] = g.getNormal(iface, nvert) * light_dir;
		return world_coord;
	}
	virtual Colorf fragmentShader(Vertex uv, Vertex barycentric) {
		Vertex vi = Vertex(varying_intensity[0], varying_intensity[1], varying_intensity[2]);
		float intensity = vi * barycentric;
		if (intensity > 5) intensity = 1;
		else if (intensity > 2) intensity = .80;
		else if (intensity > 1) intensity = .60;
		else if (intensity > .3) intensity = .45;
		else if (intensity > .15) intensity = .30;
		else intensity = 0;
		Colorf baseColor = Colorf(g.colors[varying_vert].r / 255.0, g.colors[varying_vert].g / 255.0, g.colors[varying_vert].b / 255.0);

		Colorf result = baseColor * intensity;
		return result;
	}
	virtual void setG(Geometry _g) {
		g = _g;
	}
};
struct PhongShader : public IShader {
	int varying_iface;
	int varying_vert;
	Vertex light_dir;
	Geometry g;

	virtual Vertex vertexShader(int iface, int nvert) {
		varying_iface = iface;
		varying_vert = nvert;
		Vertex world_coord = g.verts[g.faces[iface][nvert] - 1];
		return world_coord;
	}
	virtual Colorf fragmentShader(Vertex uv, Vertex barycentric) {
		Vertex N =
			g.getNormal(varying_iface, 0) * barycentric[0] +
			g.getNormal(varying_iface, 1) * barycentric[1] +
			g.getNormal(varying_iface, 2) * barycentric[2];
		N.normalize();
		Vertex l = light_dir.normalize();
		Vertex r = (N * (N * l * 2.0f) - l).normalize();	// Reflected Vector
		float spec = pow(max(r.z, 0.0f), 100);
		float diffuse = max(0.0f, N * l);
		//Colorf base(0.1, 0.4, 0.7);
		Colorf base(0, 0, 1);
		Colorf baseColor = Colorf(g.colors[varying_vert].r / 255.0, g.colors[varying_vert].g / 255.0, g.colors[varying_vert].b / 255.0);
		Colorf ambient(0.02, 0.1, 0.4);
		Colorf result;
		result.r = min<float>(ambient.r + baseColor.r * (diffuse + 2 * spec), 1);
		result.g = min<float>(ambient.g + baseColor.g * (diffuse + 2 * spec), 1);
		result.b = min<float>(ambient.b + baseColor.b * (diffuse + 2 * spec), 1);

		return result;
	}
	virtual void setG(Geometry _g) {
		g = _g;
	}
};
struct HalfToneShader : public IShader {
	float varying_intensity[3];
	int varying_vert;
	Vertex light_dir;
	Geometry g;

	virtual Vertex vertexShader(int iface, int nvert) {
		varying_vert = nvert;
		Vertex world_coord = g.verts[g.faces[iface][nvert] - 1];
		varying_intensity[nvert] = g.getNormal(iface, nvert) * light_dir;
		return world_coord;
	}
	virtual Colorf fragmentShader(Vertex uv, Vertex barycentric) {
		Vertex vi = Vertex(varying_intensity[0], varying_intensity[1], varying_intensity[2]);
		float intensity = vi * barycentric;

		Vertex newuv = Vertex(uv.x / width, uv.y / height);
		Colorf result(0, 1, 0);
		newuv = newuv * 160;
		newuv.x = frac(newuv.x);
		newuv.y = frac(newuv.y);
		newuv.z = frac(newuv.z);

		result = Colorf(newuv.x, newuv.y, 0.0);
		result = Colorf(circle(newuv, .9 * smoothstep(0.01, 10, intensity)));

		return result;
	}
	virtual void setG(Geometry _g) {
		g = _g;
	}
};






// -------------- Setup Functions -------------- 
void setupGlut() {
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0, 0, 0, 1);
	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	gluOrtho2D(0, width, 0, height);
}
void setupFile(const char* filename) {
	ifstream in(filename);
	//in.open(filename, ifstream::in);
	if (in.fail())
		return;
	string line;

	getline(in, line); // Get Geometry Nums
	int geoNums = stoi(line);
	int geoIndex = geoNums;
	cout << geoNums << endl;

	getline(in, line); // skip a blank line

	getline(in, line);
	int vertsNum = stoi(line);
	int colorNum = vertsNum;
	cout << vertsNum << endl;
	int faceNum = -1;
	int specsNum = -1;
	while (!in.eof() && geoIndex > 0) {
		if (vertsNum > 0) {
			//cout << "read vert" << vertsNum << endl;
			getline(in, line);
			vector<string> result;
			string_split(line, ' ', result);
			verts.push_back(Vertex(stof(result[0]), stof(result[1]), stof(result[2])));
			vertsNum--;
		}
		else if (colorNum > 0) {
			//cout << "read color" << colorNum << endl;
			getline(in, line);
			vector<string> result;
			string_split(line, ' ', result);
			colors.push_back(Color(stoi(result[0]), stoi(result[1]), stoi(result[2])));
			colorNum--;
		}
		else if (colorNum == 0 && faceNum == -1 && specsNum == -1) {
			getline(in, line);
			faceNum = stoi(line);
			specsNum = stoi(line);
		}
		else if (faceNum > 0) {
			//cout << "Read Face: " << faceNum << endl;
			getline(in, line);
			vector<string> result;
			string_split(line, ' ', result);
			faces.push_back(Face(stoi(result[0]), stoi(result[1]), stoi(result[2])));
			faceNum--;
		}
		else if (specsNum > 0) {
			//cout << "Read Spec: " << specsNum << endl;
			getline(in, line);
			vector<string> result;
			string_split(line, ' ', result);
			specs.push_back(stoi(result[0]));
			specsNum--;
		}
		else {
			geoIndex--;
			if (geoIndex > 0) {
				getline(in, line);	// skip a blank line
				getline(in, line);
				vertsNum = stoi(line);
				colorNum = vertsNum;
				faceNum = -1;
				specsNum = -1;
				//geos.push_back(Geometry(verts, colors, faces, specs));
				verts.clear();
				colors.clear();
				faces.clear();
				specs.clear();
			}
		}
	}

	cout << "Face Num: " << faces.size() << " Specs Num: " << specs.size() << endl;
}


// -------------- Draw Primitives Functions (Deprecated) -------------- 
void draw_pixel(int x, int y) {	// Deprecated (Use frameBufferR/G/B instead)
	glBegin(GL_POINTS);
	glVertex2i(x, y);	// Basic of EVERYTHING!!!!!
	glEnd();
}

void draw_line_bresenham(int p1x, int p1y, int p2x, int p2y) {
	int dx, dy, i, e;
	int incx, incy, inc1, inc2;
	int x, y;

	dx = p2x - p1x;
	dy = p2y - p1y;

	if (dx < 0) dx = -dx;
	if (dy < 0) dy = -dy;
	incx = 1;
	if (p2x < p1x) incx = -1;
	incy = 1;
	if (p2y < p1y) incy = -1;
	x = p1x; y = p1y;
	if (dx > dy) {
		draw_pixel(x, y);
		e = 2 * dy - dx;
		inc1 = 2 * (dy - dx);
		inc2 = 2 * dy;
		for (i = 0; i < dx; i++) {
			if (e >= 0) {
				y += incy;
				e += inc1;
			}
			else
				e += inc2;
			x += incx;
			draw_pixel(x, y);
		}

	}
	else {
		draw_pixel(x, y);
		e = 2 * dx - dy;
		inc1 = 2 * (dx - dy);
		inc2 = 2 * dx;
		for (i = 0; i < dy; i++) {
			if (e >= 0) {
				x += incx;
				e += inc1;
			}
			else
				e += inc2;
			y += incy;
			draw_pixel(x, y);
		}
	}
}
void draw_line(int p1x, int p1y, int p2x, int p2y) {
	bool steep = false;
	if (abs(p2x - p1x) < abs(p2y - p1y)) {
		swap(p1x, p1y);
		swap(p2x, p2y);
		steep = true;
	}
	if (p1x > p2x) {
		swap(p1x, p2x);
		swap(p1y, p2y);
	}
	for (float x = p1x; x < p2x; x++) {
		float t = (x - p1x) / (float)(p2x - p1x);
		int y = lerp(p1y, p2y, t);
		if (steep) {
			draw_pixel(y, x);
		}
		else {
			draw_pixel(x, y);
		}
	}
}
void draw_triangle_noFill(Vertex v1, Vertex v2, Vertex v3) {
	draw_line(v1.x, v1.y, v2.x, v2.y);
	draw_line(v2.x, v2.y, v3.x, v3.y);
	draw_line(v1.x, v1.y, v3.x, v3.y);
}
void draw_triangle_old(Vertex v1, Vertex v2, Vertex v3) {	// Line sweeping method
	if (v1.y > v2.y)
		swap(v1, v2);
	if (v1.y > v3.y)
		swap(v1, v3);
	if (v2.y > v3.y)
		swap(v2, v3);
	int total_height = v3.y - v1.y;
	for (int y = v1.y; y <= v2.y; y++) {
		int lower_segment_height = v2.y - v1.y + 1;
		float deltaA = (float)(y - v1.y) / total_height;
		float deltaB = (float)(y - v1.y) / lower_segment_height;
		Vertex A = v1 + (v3 - v1) * deltaA;
		Vertex B = v1 + (v2 - v1) * deltaB;
		if (A.x > B.x)
			swap(A, B);
		for (int x = A.x; x <= B.x; x++) {
			draw_pixel(x, y);
		}
	}
	for (int y = v2.y; y <= v3.y; y++) {
		int upper_segment_height = v3.y - v2.y + 1;
		float deltaA = (float)(y - v1.y) / total_height;
		float deltaB = (float)(y - v2.y) / upper_segment_height;
		Vertex A = v1 + (v3 - v1) * deltaA;
		Vertex B = v2 + (v3 - v2) * deltaB;
		if (A.x > B.x)
			swap(A, B);
		for (int x = A.x; x <= B.x; x++) {
			draw_pixel(x, y);
		}
	}
}
void draw_triangle(Vertex t0, Vertex t1, Vertex t2) {
	if (t0.y == t1.y && t0.y == t2.y) return;
	if (t0.y > t1.y) swap(t0, t1);
	if (t0.y > t2.y) swap(t0, t2);
	if (t1.y > t2.y) swap(t1, t2);
	int total_height = t2.y - t0.y;
	for (int i = 0; i < total_height; i++) {
		bool second_half = i > t1.y - t0.y || t1.y == t0.y;
		int segment_height = second_half ? t2.y - t1.y : t1.y - t0.y;
		float alpha = (float)i / total_height;
		float beta = (float)(i - (second_half ? t1.y - t0.y : 0)) / segment_height;
		Vertex A = t0 + (t2 - t0) * alpha;
		Vertex B = second_half ? t1 + (t2 - t1) * beta : t0 + (t1 - t0) * beta;
		if (A.x > B.x) swap(A, B);
		for (int j = A.x; j <= B.x; j++) {
			draw_pixel(j, t0.y + i);
		}
	}
}
void draw_triangle_intensity(Vertex t0, Vertex t1, Vertex t2, float ity0, float ity1, float ity2) {
	if (t0.y == t1.y && t0.y == t2.y) return;
	if (t0.y > t1.y) { swap(t0, t1); swap(ity0, ity1); }
	if (t0.y > t2.y) { swap(t0, t2); swap(ity0, ity2); }
	if (t1.y > t2.y) { swap(t1, t2); swap(ity1, ity2); }

	int total_height = t2.y - t0.y;
	for (int i = 0; i < total_height; i++) {
		bool second_half = i > t1.y - t0.y || t1.y == t0.y;
		int segment_height = second_half ? t2.y - t1.y : t1.y - t0.y;
		float alpha = (float)i / total_height;
		float beta = (float)(i - (second_half ? t1.y - t0.y : 0)) / segment_height; // be careful: with above conditions no division by zero here
		Vertex A = t0 + Vertex(t2 - t0) * alpha;
		Vertex B = second_half ? t1 + Vertex(t2 - t1) * beta : t0 + Vertex(t1 - t0) * beta;
		float ityA = ity0 + (ity2 - ity0) * alpha;
		float ityB = second_half ? ity1 + (ity2 - ity1) * beta : ity0 + (ity1 - ity0) * beta;
		if (A.x > B.x) {
			swap(A, B);
			swap(ityA, ityB);
		}
		for (int j = A.x; j <= B.x; j++) {
			float phi = B.x == A.x ? 1. : (float)(j - A.x) / (B.x - A.x);
			Vertex P = Vertex(A) + Vertex(B - A) * phi;
			float ityP = ityA + (ityB - ityA) * phi;
			int idx = P.x + P.y * width;
			if (P.x >= width || P.y >= height || P.x < 0 || P.y < 0) continue;

			if (ityP > 4) ityP = 1;
			else if (ityP > 0.7) ityP = 0.7;
			else if (ityP > 0.3) ityP = 0.4;
			else if (ityP > 0.1) ityP = 0.2;
			else ityP = 0;
			//ityP *= 0.15;
			glColor3f(ityP, ityP, ityP);
			draw_pixel(P.x, P.y);

		}
	}
}


// -------------- Draw Primitives Functions (Currently Using) -------------- 

Vertex barycentric(Vertex A, Vertex B, Vertex C, Vertex P) {
	Vertex s[2];
	for (int i = 2; i--; ) {
		s[i].x = C[i] - A[i];
		s[i].y = B[i] - A[i];
		s[i].z = A[i] - P[i];
	}
	Vertex u = crossProduct(s[0], s[1]);
	if (abs(u[2]) > 1e-2)
		return Vertex(1.f - (u.x + u.y) / u.z, u.y / u.z, u.x / u.z);
	return Vertex(-1, 1, 1);
}
void draw_triangle_barycentric(Vertex* pts, IShader& shader) {
	Vertex bboxmin(numeric_limits<int>::max(), numeric_limits<int>::max());
	Vertex bboxmax(-numeric_limits<int>::max(), -numeric_limits<int>::max());
	Vertex clamp(width - 1, height - 1);
	for (int i = 0; i < 3; i++) {
		bboxmin.x = max(0.0f, min(bboxmin.x, pts[i].x));
		bboxmax.x = min(clamp.x, max(bboxmax.x, pts[i].x));
		bboxmin.y = max(0.0f, min(bboxmin.y, pts[i].y));
		bboxmax.y = min(clamp.y, max(bboxmax.y, pts[i].y));
	}
	Vertex P;
	for (P.x = bboxmin.x; P.x <= bboxmax.x; P.x++) {
		for (P.y = bboxmin.y; P.y <= bboxmax.y; P.y++) {
			Vertex bc_screen = barycentric(pts[0], pts[1], pts[2], P);
			if (bc_screen.x < 0 || bc_screen.y < 0 || bc_screen.z < 0) continue;
			P.z = 0;
			for (int i = 0; i < 3; i++)
				P.z += pts[i].z * bc_screen[i];

			Colorf col = shader.fragmentShader(P, bc_screen);

			frameBufferR[int(P.x)][int(P.y)] = col.r * 255;
			frameBufferG[int(P.x)][int(P.y)] = col.g * 255;
			frameBufferB[int(P.x)][int(P.y)] = col.b * 255;
			//glColor3f(col.r, col.g, col.b);
			//draw_pixel(P.x, P.y);
		}
	}
}


class Zone {
public:
	Model* model;
	vector<vector<Face>> zone_faces;
	vector<Geometry> zone_g;

	int axis;
	float gScal = 1.0;
	float gOffsetX = width / 2;
	float gOffsetY = height / 2;

	Zone(Model* _model, int _axis) {
		model = _model;
		axis = _axis;


		for (int i = 0; i < model->geos.size(); i++) {
			model->geos[i].sortFaces(axis);
		}
		zone_g.assign(model->geos.begin(), model->geos.end());
	}

	Vertex world2screen(Vertex v, float _gScal, float _gOffsetX, float _gOffsetY) {
		if (axis == 0)	// xy - z
			return Vertex(int((v.x - 0.5) * width * _gScal + _gOffsetX), int((v.y - 0.5) * height * _gScal + _gOffsetY), int(v.z));
		if (axis == 1)	// xz - y
			return Vertex(int((v.x - 0.5) * width * _gScal + _gOffsetX), int((v.z - 0.5) * height * _gScal + _gOffsetY), int(v.y));
		if (axis == 2)	// zy - x
			return Vertex(int((v.z - 0.5) * width * _gScal + _gOffsetX), int((v.y - 0.5) * height * _gScal + _gOffsetY), int(v.x));
	}
	void drawZone(float scale, int zoneX, int zoneY, IShader& shader) {
		for (int g = 0; g < model->geos.size(); g++) {
			shader.setG(zone_g[g]);	// Did some trick
			for (int i = 0; i < model->getGeo()[g].faces.size(); i++) {
				Vertex screen_coords[3];
				for (int j = 0; j < 3; j++) {
					screen_coords[j] = world2screen(shader.vertexShader(i, j), scale, zoneX, zoneY);
				}
				draw_triangle_barycentric(screen_coords, shader);
			}
		}
	}

};

Model* model;

Zone* cube_xy;
Zone* cube_xz;
Zone* cube_yz;
Zone* bunny_xy;
Zone* bunny_xz;
Zone* bunny_yz;
Zone* ven_xy;
Zone* ven_xz;
Zone* ven_yz;


GourandShader gourand;
PhongShader phong;
HalfToneShader halftone;
ToonShader toon;



// -------------- Display Function -------------- 

void myDisplay() {

	//// ----------- Draw WireFrame -------------------
	//glColor3f(1, 0, 1);
	//for (int i = 0; i < faces.size(); i++) {
	//	for (int j = 0; j < 3; j++) {
	//		Vertex v1 = verts[faces[i][j] - 1];
	//		Vertex v2 = verts[faces[i][(j + 1) % 3] - 1];
	//		draw_line(int(v1.x * width), int(v1.y * height), int(v2.x * width), int(v2.y * height));
	//	}
	//}


	//// ----------- Draw Colored Triangles -------------
	//for (int i = 0; i < faces.size(); i++) {
	//	Vertex screen_coords[3];
	//	for (int j = 0; j < 3; j++) {
	//		screen_coords[j] = world2screen(verts[faces[i][j] - 1]);
	//	}
	//	glColor3f((rand() % 100) / 100.0, (rand() % 100) / 100.0, (rand() % 100) / 100.0);
	//	draw_triangle(screen_coords[0], screen_coords[1], screen_coords[2]);
	//}



	//// ---------- Draw Lighted Triangles ---------------
	//Vertex light_dir(sin(timer), 0, 1);
	//for (int i = 0; i < faces.size(); i++) {
	//	Vertex screen_coords[3];
	//	Vertex world_coords[3];
	//	for (int j = 0; j < 3; j++) {
	//		screen_coords[j] = world2screen(verts[faces[i][j] - 1]);
	//		world_coords[j] = verts[faces[i][j] - 1];
	//	}
	//	Vertex N = crossProduct((world_coords[1] - world_coords[0]), (world_coords[2] - world_coords[0]));
	//	N.normalize();
	//	float intensity = N * light_dir;
	//	if (intensity > 0) {
	//		glColor3f(intensity, intensity, intensity);
	//		draw_triangle(screen_coords[0], screen_coords[1], screen_coords[2]);
	//	}
	//	else {
	//		glColor3f(0, 0, 0);
	//		draw_triangle(screen_coords[0], screen_coords[1], screen_coords[2]);
	//	}
	//}



	//// -------- Draw Barycentric Triangles ---------------
	//Vertex light_dir(sin(timer), 0, 1);
	//for (int i = 0; i < faces.size(); i++) {
	//	Vertex screen_coords[3];
	//	Vertex world_coords[3];

	//	for (int j = 0; j < 3; j++) {
	//		screen_coords[j] = world2screen(verts[faces[i][j] - 1]);
	//		world_coords[j] = verts[faces[i][j] - 1];
	//	}

	//	Vertex N = crossProduct((world_coords[1] - world_coords[0]), (world_coords[2] - world_coords[0]));
	//	N.normalize();
	//	float intensity = N * light_dir;
	//	if (intensity > 0) {
	//		glColor3f(intensity, intensity, intensity);
	//		draw_triangle_barycentric(screen_coords);
	//	}
	//	else {
	//		glColor3f(0, 0, 0);
	//		draw_triangle_barycentric(screen_coords);
	//	}
	//}



	//// ----------- Gouraud shading --------------

	//Vertex light_dir(sin(timer) * 3, 0, 1);

	//for (int i = 0; i < faces.size(); i++) {
	//	Vertex screen_coords[3];
	//	float intensity[3];
	//	for (int j = 0; j < 3; j++) {
	//		screen_coords[j] = world2screen(verts[faces[i][j] - 1]);
	//		intensity[j] = norms[getNormIndex(i, j) - 1] * light_dir;
	//	}
	//	draw_triangle_intensity(screen_coords[0], screen_coords[1], screen_coords[2], intensity[0], intensity[1], intensity[2]);
	//}





	//// ------------ Universal Shading ---------------
	//Vertex light_dir(sin(timer), 0, 1);

	//GourandShader gourandShader;
	//gourandShader.light_dir = light_dir;
	//ToonShader toonShader;
	//toonShader.light_dir = light_dir;
	//PhongShader phongShader;
	//phongShader.light_dir = light_dir;

	//for (int i = 0; i < faces.size(); i++) {
	//	Vertex screen_coords[3];
	//	for (int j = 0; j < 3; j++) {
	//		screen_coords[j] = phongShader.vertexShader(i, j);
	//	}
	//	draw_triangle_barycentric(screen_coords, phongShader);
	//}



	//// ------------ Multiple Model Loading ---------------
	//Vertex light_dir(sin(timer) * 1.5, cos(timer), 1);

	//GourandShader gourandShader;
	//gourandShader.light_dir = light_dir;
	//ToonShader toonShader;
	//toonShader.light_dir = light_dir;
	//PhongShader phongShader;
	//phongShader.light_dir = light_dir;
	//HalfToneShader halftoneShader;
	//halftoneShader.light_dir = light_dir;

	//for (int g = 0; g < model->geos.size(); g++) {
	//	halftoneShader.g = model->getGeo()[g];
	//	for (int i = 0; i < model->getGeo()[g].faces.size(); i++) {
	//		Vertex screen_coords[3];
	//		for (int j = 0; j < 3; j++) {
	//			screen_coords[j] = world2screen(halftoneShader.vertexShader(i, j));
	//		}
	//		draw_triangle_barycentric(screen_coords, halftoneShader);
	//	}
	//}

	// --------------- Draw Scene --------------------------

	if (!hasDrawn) {
		Vertex light_dir(1.5, 1, 1);

		if (modelActiveBtn == 0) {
			light_dir = Vertex(0, 0.5, 1);
		}
		else if (modelActiveBtn == 1) {
			light_dir = Vertex(1.5, 1, 1);
		}
		gourand.light_dir = light_dir;
		phong.light_dir = light_dir;
		halftone.light_dir = light_dir;
		toon.light_dir = light_dir;


		if (shaderActiveBtn == 0) {
			if (modelActiveBtn == 0) {
				bunny_xy->drawZone(2, width / 4, height / 4 * 2, phong);
				bunny_xz->drawZone(2, width / 4 * 3, height / 4 * 3, phong);
				bunny_yz->drawZone(2, width / 4, 0, phong);
			}
			else if (modelActiveBtn == 1) {
				cube_xy->drawZone(0.9, width / 7, height / 10 * 8.5, phong);
				cube_xz->drawZone(0.9, width / 8 * 5, height / 10 * 8.5, phong);
				cube_yz->drawZone(0.9, width / 7, height / 10 * 3.5, phong);
			}
			else if (modelActiveBtn == 2) {
				ven_xy->drawZone(1, width / 4, height / 4 * 3, phong);
				ven_xz->drawZone(1, width / 4 * 3, height / 4 * 3, phong);
				ven_yz->drawZone(1, width / 4, height / 4, phong);
			}
		}
		else if (shaderActiveBtn == 1) {
			if (modelActiveBtn == 0) {
				bunny_xy->drawZone(2, width / 4, height / 4 * 2, gourand);
				bunny_xz->drawZone(2, width / 4 * 3, height / 4 * 3, gourand);
				bunny_yz->drawZone(2, width / 4, 0, gourand);
			}
			else if (modelActiveBtn == 1) {
				cube_xy->drawZone(0.9, width / 7, height / 10 * 8.5, gourand);
				cube_xz->drawZone(0.9, width / 8 * 5, height / 10 * 8.5, gourand);
				cube_yz->drawZone(0.9, width / 7, height / 10 * 3.5, gourand);
			}
			else if (modelActiveBtn == 2) {
				ven_xy->drawZone(1, width / 4, height / 4 * 3, gourand);
				ven_xz->drawZone(1, width / 4 * 3, height / 4 * 3, gourand);
				ven_yz->drawZone(1, width / 4, height / 4, gourand);
			}
		}
		else if (shaderActiveBtn == 2) {
			if (modelActiveBtn == 0) {
				bunny_xy->drawZone(2, width / 4, height / 4 * 2, halftone);
				bunny_xz->drawZone(2, width / 4 * 3, height / 4 * 3, halftone);
				bunny_yz->drawZone(2, width / 4, 0, halftone);
			}
			else if (modelActiveBtn == 1) {
				cube_xy->drawZone(0.9, width / 7, height / 10 * 8.5, halftone);
				cube_xz->drawZone(0.9, width / 8 * 5, height / 10 * 8.5, halftone);
				cube_yz->drawZone(0.9, width / 7, height / 10 * 3.5, halftone);
			}
			else if (modelActiveBtn == 2) {
				ven_xy->drawZone(1, width / 4, height / 4 * 3, halftone);
				ven_xz->drawZone(1, width / 4 * 3, height / 4 * 3, halftone);
				ven_yz->drawZone(1, width / 4, height / 4, halftone);
			}
		}
		else if (shaderActiveBtn == 3) {
			if (modelActiveBtn == 0) {
				bunny_xy->drawZone(2, width / 4, height / 4 * 2, toon);
				bunny_xz->drawZone(2, width / 4 * 3, height / 4 * 3, toon);
				bunny_yz->drawZone(2, width / 4, 0, toon);
			}
			else if (modelActiveBtn == 1) {
				cube_xy->drawZone(0.9, width / 7, height / 10 * 8.5, toon);
				cube_xz->drawZone(0.9, width / 8 * 5, height / 10 * 8.5, toon);
				cube_yz->drawZone(0.9, width / 7, height / 10 * 3.5, toon);
			}
			else if (modelActiveBtn == 2) {
				ven_xy->drawZone(1, width / 4, height / 4 * 3, toon);
				ven_xz->drawZone(1, width / 4 * 3, height / 4 * 3, toon);
				ven_yz->drawZone(1, width / 4, height / 4, toon);
			}
		}

		hasDrawn = true;
	}

	// ---------------- Display the frameBuffers -------------------
	glBegin(GL_POINTS);
	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			glColor3f(frameBufferR[i][j] / 255., frameBufferG[i][j] / 255., frameBufferB[i][j] / 255.);
			glVertex2i(i, j);
		}
	}
	glEnd();


	// ---------------- Draw UI ---------------------------
	for (int i = 0; i < shader_btns.size(); i++) {
		shader_btns[i]->drawBtn();
		if (shader_btns[i]->state == PRESSED) {
			shaderActiveBtn = i;
		}
	}
	for (int i = 0; i < model_btns.size(); i++) {
		model_btns[i]->drawBtn();
		if (model_btns[i]->state == PRESSED) {
			modelActiveBtn = i;
		}
	}



	//glutSwapBuffers();
	glFlush();
	glutPostRedisplay();
}


void myTimerFunc(int) {
	glutPostRedisplay();
	glutTimerFunc(1000 / 60, myTimerFunc, 0);
	cout << sin(timer) << endl;
	timer += 0.1;
}

int main(int argc, char** argv) {
	// Setup Window
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowSize(width, height);
	glutInitWindowPosition(0, 0);
	glutCreateWindow("GlutHW");


	initUI();
	setupGlut();

	// ------------- Load Models -----------------
	model = new Model("cube.txt");
	for (int i = 0; i < model->geos.size(); i++) {
		model->geos[i].sortFaces(0);
		model->geos[i].calcNormals();
	}
	cube_xy = new Zone(model, 0);
	cube_xz = new Zone(model, 1);
	cube_yz = new Zone(model, 2);


	model = new Model("bunny.txt");
	for (int i = 0; i < model->geos.size(); i++) {
		model->geos[i].sortFaces(0);
		model->geos[i].calcNormals();
	}
	bunny_xy = new Zone(model, 0);
	bunny_xz = new Zone(model, 1);
	bunny_yz = new Zone(model, 2);


	model = new Model("Beethoven.txt");
	for (int i = 0; i < model->geos.size(); i++) {
		model->geos[i].sortFaces(0);
		model->geos[i].calcNormals();
	}
	ven_xy = new Zone(model, 0);
	ven_xz = new Zone(model, 1);
	ven_yz = new Zone(model, 2);



	glutDisplayFunc(myDisplay);

	//glutTimerFunc(0, myTimerFunc, 0);

	glutMouseFunc(mouse_hit);
	glutPassiveMotionFunc(mouse_move);

	glutMainLoop();

	return 0;
}

IShader::~IShader()
{
}
