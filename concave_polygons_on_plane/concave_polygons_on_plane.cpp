#ifdef _WIN32
#define NOMINMAX
#include <Windows.h>
#define _USE_MATH_DEFINES
#define copysign _copysign
double time() {
    LARGE_INTEGER  c, f;
    QueryPerformanceFrequency(&f);
    QueryPerformanceCounter(&c);
    return (double)c.QuadPart / (double)f.QuadPart;
}
#else
#include <time.h>
double time() {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return 1.0 * ts.tv_sec + ts.tv_nsec / 1000000000.0;
}
#endif

#include <cassert>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <GL/glut.h>

struct Vector {
    float x, y;
    Vector(float X, float Y): x(X), y(Y) {}
    Vector(){}
    float length() const { return sqrtf(length2()); }
    float length2() const { return x*x + y*y;  }
    Vector operator*(float a) const { return Vector(x*a, y*a); }
    Vector operator-(const Vector& b) const { return Vector(x - b.x, y - b.y); }
    Vector operator+(const Vector& b) const { return Vector(x + b.x, y + b.y); }
};
typedef Vector Point;


float cross(const Vector&a, const Vector&b) { return a.x * b.y - a.y * b.x; }
float dot(const Vector&a, const Vector&b) { return a.x * b.x + a.y * b.y; }

struct Color {
    float r, g, b, a;
    Color(float R, float G, float B, float A = 1): r(R), g(G), b(B), a(A) {}
    Color() {}
};

typedef std::vector<Point> PointArray;
typedef std::vector<Color> ColorArray;

struct Stipple {
    PointArray Vertices;
    ColorArray Colors;
    Stipple(const PointArray& contour) {
        BuildStipple(-0.075, contour);
        SetColor(Color(0.5, 0.5, 0.5));
    }

    void SetColor(Color nearColor, Color farColor) {
        Colors.resize(Vertices.size());
        for (size_t i = 0; i < Colors.size(); ++i)
            Colors[i] = i % 2 ? nearColor : farColor;
    }
    void SetColor(Color basicColor) {
        SetColor(Color(basicColor.r, basicColor.g, basicColor.b, 1), Color(basicColor.r, basicColor.g, basicColor.b, 0.125));
    }

    void BuildStipple(float distance, const PointArray& a) {
        size_t sz = a.size();
        PointArray& ll = Vertices;
        ll.clear();
        for (size_t i = 0; i < sz; ++i) {
            size_t prev = (i + sz - 1) % sz;
            size_t next = (i + 1) % sz;
            Vector v0 = a[i] - a[prev];
            Vector v1 = a[next] - a[i];
            float c = cross(v0, v1);
            Vector r0 = Vector(v0.y, -v0.x) * copysign(1 / v0.length(), distance);
            Vector r1 = Vector(v1.y, -v1.x) * copysign(1 / v1.length(), distance);
            float factor = fabsf(distance);
            Point p0 = a[i] + r0 * factor;
            Point p1 = a[i] + r1 * factor;
            if (c * distance > 0) {
                Vector d = p0 + (p1 - p0) * 0.5 - a[i];
                Vector p = a[i] + d * (factor / d.length());
                ll.push_back(p0);
                ll.push_back(a[i]);
                ll.push_back(p);
                ll.push_back(a[i]);
                ll.push_back(p1);
                ll.push_back(a[i]);

            } else {
                Vector d = p0 + (p1 - p0) * 0.5 - a[i];
                Vector p = a[i] + d * (factor / dot(d, r0));
                ll.push_back(p);
                ll.push_back(a[i]);
            }
        }
        ll.push_back(ll[0]);
        ll.push_back(ll[1]);
    }
    void Draw() const {
        glEnable( GL_POLYGON_STIPPLE);
        glVertexPointer(2, GL_FLOAT, sizeof(Vertices[0]), &Vertices[0].x);
        glColorPointer(4, GL_FLOAT, sizeof(Colors[0]), &Colors[0].r);
        glEnableClientState(GL_VERTEX_ARRAY);
        glEnableClientState(GL_COLOR_ARRAY);
        glDrawArrays(GL_QUAD_STRIP, 0, (GLsizei)Vertices.size());
        glDisableClientState(GL_VERTEX_ARRAY);
        glDisableClientState(GL_COLOR_ARRAY);
        glDisable(GL_POLYGON_STIPPLE);
    }
};

struct TransformMatrix { GLfloat m[16]; };
const TransformMatrix IdentityMatrix = {{1, 0, 0, 0,  0, 1, 0, 0,  0, 0, 1, 0,  0, 0, 0, 1}};

PointArray ReverseIf(PointArray a, bool reverse) { if (reverse) std::reverse(a.begin(), a.end()); return a; }

struct Figure {
    bool Inverse;
    PointArray Vertices;
    Stipple InnerStipple;
    Stipple SecondInnerStipple;
    TransformMatrix Transform;
    bool Active;

    Figure(const PointArray& p, bool inv, Color c, bool active = true)
        : Inverse(inv)
        , Vertices(ReverseIf(p, Inverse))
        , InnerStipple(Vertices)
        , SecondInnerStipple(InnerStipple)
        , Transform(IdentityMatrix)
        , Active(active) {
          InnerStipple.SetColor(c);
          SecondInnerStipple.SetColor(Color(0.1, 0.1, 0.1));
        }
    void DrawStipple(bool second) const { (second ? SecondInnerStipple : InnerStipple).Draw(); }
};

void DrawContour(const PointArray& a, GLenum mode) {
    glVertexPointer(2, GL_FLOAT, sizeof(a[0]), &a[0].x);
    glEnableClientState(GL_VERTEX_ARRAY);
    glDrawArrays(mode, 0, (GLsizei)a.size());
    glDisableClientState(GL_VERTEX_ARRAY);
}

void DrawRectangle() {
    Point v[4] = {Point(-1, -1), Point(1, -1), Point(1, 1), Point(-1, 1)};
    glVertexPointer(2, GL_FLOAT, sizeof(Point), &v[0].x);
    glEnableClientState(GL_VERTEX_ARRAY);
    glDrawArrays(GL_QUADS, 0, 4);
    glDisableClientState(GL_VERTEX_ARRAY);
}

float lerp(float a, float b, float t) { return a + (b - a) * t; }

PointArray CreateStar(int n, int m, float inner, float outer, float offset) {
    assert(m >= 2);
    PointArray a;
    a.reserve(n * m);
    int centerIndex = (m + 1) / 2;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            float fr = (j <= centerIndex ? j : 2 * centerIndex - j) / (float)centerIndex;
            float r = lerp(inner, outer, fr);
            float angle = (2 * M_PI) / (m * n) * (j + i * m) + offset * fr;
            a.push_back(Vector(cos(angle) * r , sin(angle) * r));
        }
    }
    return a;
}

struct FigureHolder {
    std::vector<Figure*> Figures;
    void DrawStripple(size_t idx, bool colorless) const {
        glPushMatrix();
        glMultMatrixf(Figures[idx]->Transform.m);
        Figures[idx]->DrawStipple(colorless);
        glPopMatrix();
    }
    void DrawAllStripples(bool colorless = false) const {
        for (size_t i = 0; i < Figures.size(); ++i)
            if (Figures[i]->Active)
                DrawStripple(i, colorless);
    }
    void DrawContour(size_t idx, GLenum mode) const{
        glPushMatrix();
        glMultMatrixf(Figures[idx]->Transform.m);
        ::DrawContour(Figures[idx]->Vertices, mode);
        glPopMatrix();
    }
    void DrawAllContours(GLenum mode) const{
        for (size_t i = 0; i < Figures.size(); ++i)
            if (Figures[i]->Active)
                DrawContour(i, mode);
    }
};

unsigned PrepareIntersecionMask(const FigureHolder& p) {
    size_t n = p.Figures.size();
    unsigned ref = 0;
    glEnable(GL_STENCIL_TEST);
    glStencilMask(-1);
    glClear(GL_STENCIL_BUFFER_BIT);
    glStencilFunc(GL_NEVER, 0, -1);
    glStencilOp(GL_INVERT, GL_INVERT, GL_INVERT);
    for (size_t i = 0; i < n; ++i)
        if (p.Figures[i]->Active) {
            glStencilMask(1 << i);
            p.DrawContour(i, GL_TRIANGLE_FAN);
            ref |= !p.Figures[i]->Inverse << i;
        }
    glDisable(GL_STENCIL_TEST);
    return ref;
}


const int N = 6;
Figure f[N] = {
    Figure(CreateStar(5, 3, 0.4, 0.9, 0), false, Color(0, 1, 0)),
    Figure(CreateStar(2, 3, 0.275, 1.1, 0), false, Color(0.3, 0.5, 1)),
    Figure(CreateStar(4, 3, 0.05, 0.1, 0), true, Color(1, 0, 0)),
    Figure(CreateStar(11, 30, 0.5, 0.7, 0), false, Color(1, 0.7, 0), false),
    Figure(CreateStar(5, 40, 0.4, 0.9, 0.5), false, Color(0.8, 0.2, 1), false),
    Figure(CreateStar(5, 40, 0.4, 0.9, -0.7), false, Color(0.8, 1, 0.2), false),

};
int width = 1024;
int height = 512;
bool drawSourceStipple = true;
bool drawTargetStipple = true;
FigureHolder holder;

void animate() {
    double t = time();
    glPushMatrix();

    glLoadIdentity();
    glTranslatef(cos(t) * 0.25, sin(t) * 0.25, 0);
    glGetFloatv(GL_MODELVIEW_MATRIX, f[4].Transform.m);

    glLoadIdentity();
    glRotatef(t * 5, 0, 0, -1);
    glGetFloatv(GL_MODELVIEW_MATRIX, f[5].Transform.m);

    glPopMatrix();
}

void display() {
    glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-1, 1, -1, 1, -1, 1);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glViewport(0, 0, width / 2, height);
    if (drawSourceStipple)
        holder.DrawAllStripples();
    glLineWidth(1);
    glColor3f(0, 0, 0);
    holder.DrawAllContours(GL_LINE_LOOP);

    glViewport(width / 2, 0, width / 2, height);
    unsigned refMask = PrepareIntersecionMask(holder);


    glEnable(GL_STENCIL_TEST);
    glStencilFunc(GL_NOTEQUAL, refMask, -1);
    glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP);
    glColor3f(0.4, 0.4, 0.4);
    DrawRectangle();

    glStencilFunc(GL_EQUAL, refMask, -1);
    glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP);

    if (drawTargetStipple)
        holder.DrawAllStripples(true);

    glLineWidth(2);
    glColor3f(0, 0, 0);
    holder.DrawAllContours(GL_LINE_LOOP);

    glDisable(GL_STENCIL_TEST);

    glutSwapBuffers();
}

void reshape(int w, int h) {
    if (h == 0) h = 1;
    width = w;
    height = h;
    glViewport(0, 0, w, h);
}

void initialize() {
    for (int i = 0; i < N; ++i)
        holder.Figures.push_back(&f[i]);

    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    {
        GLubyte m[32][4] = {{0}};
        for (int i = 0; i < 32; i++)
            for (int j = 0; j < 32; j++)
                m[i][j >> 3] |= (((i + j) >> 1) & 1) << (j & 7);
        glPolygonStipple(m[0]);
    }
}


void idle() {
    animate();
    glutPostRedisplay();
}


void keyboad(unsigned char key, int, int) {
    if (key == 27)
        exit(0);
    if (key  >= '1' && key < '1' + N) {
        Figure& c = f[key - '1'];
        c.Active = !c.Active;
    }
    if (key == 's' || key == 'S')
        drawSourceStipple = !drawSourceStipple;
    if (key == 'd' || key == 'D')
        drawTargetStipple = !drawTargetStipple;
}

int main(int argc, char* argv[]) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_STENCIL | GLUT_MULTISAMPLE);
    glutInitWindowSize(width, height);
    glutCreateWindow("Concave polygon demo");

    printf("Hit ESC key to quit;\n");
    printf("1-%d\t - toggle figure;\n", N);
    printf("s\t - toggle source stipple;\n");
    printf("d\t - toggle target stipple;\n");

    glutIdleFunc(idle);
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutKeyboardFunc(keyboad);

    initialize();
    glutMainLoop();
    return 0;
}





