
#include<iostream>
#include <graphics.h>
#include <math.h>
#define color_A 11
#define color_B 10
#define ScreenHeigh 600
#define ScreenWidth 1000
#define N 10
#define M 4
struct Point
{
	float x;
	float y;
};
class xyz
{
public:
  double x,y,z;
};
bool NotIntersectPlaneLine;
struct Vertex
{
	Point position;
};
struct Vector
{
	float x;
	float y;
};
Vector operator-(Point lhs, Point rhs)
{
	Vector result;

	result.x = lhs.x - rhs.x;
	result.y = lhs.y - rhs.y;

	return result;
}
float Edge(Point v0, Point v1, Point p)
{
	Vector a = p - v0;
	Vector b = v1 - v0;
	return a.x * b.y - a.y * b.x;
}
const int DX = 10, DY = 10,DZ=10;
const double S1 = 1.1;
const double S2 = 0.95;
const double ALPHA = 0.087;
double Zbuffer[ScreenWidth][ScreenHeigh]={};
void draw_shdow(double fig1[N][M],double fig2[N][M]);
using namespace std;
void DrawFloor(double fig[N][M],double fig2[N][M], double light[4]);
void mult(double fig[N][M], double mass[M][M]);
void move(double fig[N][M], double dx, double dy,double dz);
void scale(double fig[N][M], double S);
void rotateZ(double fig[N][M], double angle);
void rotateY(double fig[N][M], double angle);
void rotateX(double fig[N][M], double angle);
void draw(double fig1[N][M],double fig2[N][M]);
double aver(double range[N][M], int a);
void bresenhamline(int x0,int y0,int x1,int y1, int color);
void nullZbuf();
void trianglePaint(int x1,int y1,int z1,int x2,int y2,int z2, int x3,int y3,int z3, int color);
xyz CreateVector(xyz A, xyz B);
xyz VectorProduct(xyz A, xyz B);
double DotProduct(xyz A, xyz B);
void Normalize(xyz &A);
xyz PlaneIntersectLine(xyz A, xyz B, xyz C, xyz X, xyz Y);
void resetVewport();
int main() {
    int chosenfigure = 1;
	setlocale(LC_ALL, "Russian");
	initwindow(ScreenWidth, ScreenHeigh);
	double figure[N][M] = {{400,10 ,500 , 1},
						  {420, 0, 500, 1},
						  {410, 30, 500, 1},
						  {440, 10, 500, 1},
						  {430, 30, 500, 1},
						  {400,10 ,550 , 1},
						  {420, 0, 550, 1},
						  {410, 30, 550, 1},
						  {440, 10, 550, 1},
						  {430, 30, 550, 1}};
    double figure2 [N][M] = {{620,0 ,500 , 1},
						  {640, 40, 500, 1},
						  {600, 40, 500, 1},
						  {620, 0, 540, 1},
						  {640, 40, 540, 1},
						  {600,40 ,540 , 1},
						  {620, 20, 520, 1},
						  {620, 20, 520, 1},
						  {620, 20, 520, 1},
						  {620, 20, 520, 1}};
    double Lightsouce[4]={500,-1000,500,0};
	cout<< "Инструкция:\n\t\tКнопки:\t\nПеремещение\t| W - вверх, S - вниз, А - влево, D - вправо\t|\n|F - вперёд, B - назад, 1, 2 - переключение между фигурами\t|\nМасштабирование\t| Q - увеличить, E - уменьшить\t\t\t|\nПоворот вокруг\t|\t\t\t\t\t\t|\nсвоей оси\t| Z - по Z, X - пo X, C - по Y\t\t\t|\n"<<endl;
cout<<"-----------------------------------------\nВыход - esc."<<endl;
DrawFloor(figure, figure2,Lightsouce);
draw(figure,figure2);

do {

if (GetAsyncKeyState(VK_ESCAPE) & 0x8000) break;
if (GetAsyncKeyState((unsigned short)'1') & 0x8000)
    chosenfigure = 1;
if (GetAsyncKeyState((unsigned short)'2') & 0x8000)
    chosenfigure = 2;
if (GetAsyncKeyState((unsigned short)'W') & 0x8000){
    if (chosenfigure ==1)
        move(figure, 0, -DY,0);
    if (chosenfigure ==2)
        move(figure2, 0, -DY,0);
    resetVewport(),DrawFloor(figure,figure2,Lightsouce), draw(figure,figure2);
}

if (GetAsyncKeyState((unsigned short)'S') & 0x8000){
    if (chosenfigure ==1)
        move(figure, 0, DY,0);
    if (chosenfigure ==2)
       move(figure2, 0, DY,0);
    resetVewport(),DrawFloor(figure,figure2,Lightsouce), draw(figure,figure2);
}

if (GetAsyncKeyState((unsigned short)'A') & 0x8000){
    if (chosenfigure ==1)
        move(figure, -DX, 0,0);
    if (chosenfigure ==2)
        move(figure2, -DX, 0,0);
    resetVewport(),DrawFloor(figure,figure2,Lightsouce), draw(figure,figure2);
}

if (GetAsyncKeyState((unsigned short)'D') & 0x8000){
    if (chosenfigure ==1)
       move(figure, DX, 0,0);
    if (chosenfigure ==2)
        move(figure2, DX, 0,0);
    resetVewport(),DrawFloor(figure,figure2,Lightsouce), draw(figure,figure2);
}

if (GetAsyncKeyState((unsigned short)'B') & 0x8000){
    if (chosenfigure ==1)
        move(figure, 0, 0,DZ);
    if (chosenfigure ==2)
        move(figure2, 0, 0,DZ);
    resetVewport(),DrawFloor(figure,figure2,Lightsouce), draw(figure,figure2);
}

if (GetAsyncKeyState((unsigned short)'F') & 0x8000)
{
    if (chosenfigure ==1)
       move(figure, 0, 0,-DZ);
    if (chosenfigure ==2)
       move(figure2, 0, 0,-DZ);
    resetVewport(),DrawFloor(figure,figure2,Lightsouce), draw(figure,figure2);
}

if (GetAsyncKeyState((unsigned short)'Q') & 0x8000){
    if (chosenfigure ==1)
        scale(figure, S1);
    if (chosenfigure ==2)
        scale(figure2, S1);
    resetVewport(),DrawFloor(figure,figure2,Lightsouce), draw(figure,figure2);
}

if (GetAsyncKeyState((unsigned short)'E') & 0x8000){
    if (chosenfigure ==1)
        scale(figure, S2);
    if (chosenfigure ==2)
        scale(figure2, S2);
    resetVewport(),DrawFloor(figure,figure2,Lightsouce), draw(figure,figure2);
}

if (GetAsyncKeyState((unsigned short)'Z') & 0x8000){
    if (chosenfigure ==1)
        rotateZ(figure, ALPHA);
    if (chosenfigure ==2)
        rotateZ(figure2, ALPHA);
    resetVewport(),DrawFloor(figure,figure2,Lightsouce), draw(figure,figure2);
}

if (GetAsyncKeyState((unsigned short)'X') & 0x8000){
    if (chosenfigure ==1)
        rotateX(figure, ALPHA);
    if (chosenfigure ==2)
       rotateX(figure2, ALPHA);
    resetVewport(),DrawFloor(figure,figure2,Lightsouce), draw(figure,figure2);
}

if (GetAsyncKeyState((unsigned short)'C') & 0x8000){
    if (chosenfigure ==1)
        rotateY(figure, ALPHA);
    if (chosenfigure ==2)
        rotateY(figure2, ALPHA);
    resetVewport(),DrawFloor(figure,figure2,Lightsouce), draw(figure,figure2);
}
delay(10);
}while (true);
	closegraph();
	return 0;
}
void mult(double fig[N][M], double mass[M][M]) {
	double res[N][M] = { 0,0,0,0 };
	for (int k = 0; k < N; k++) {
		for (int i = 0; i< M; i++) {
			for (int j = 0; j < M; j++)
				res[k][i] += fig[k][j] * mass[j][i];
		}
	}
	for (int k = 0; k < N; k++) {
		for (int i = 0; i< M; i++)
			fig[k][i] = res[k][i];
	}
}
void move(double fig[N][M], double dx, double dy,double dz) {
	double DX_DY[M][M] = { {1, 0, 0, 0} ,
							{0, 1, 0, 0},
							{0, 0, 1, 0},
                            {dx, dy, dz,1} };
	mult(fig, DX_DY);
}
void DrawFloor(double fig[N][M],double fig2[N][M], double light[4]){
    double Floor[4][4]={{0,600,0,0},{1000,600,0,0},{1000,450,1000,0},{0,450,1000,0}};
    double vect[4] = {0,0,0,0};
    double shadow[N][M] = {0};
    double shadow2[N][M] = {0};
    xyz pt;
    xyz A;
    xyz B;
    xyz C;
    xyz X;
    xyz Y;
    A.x=Floor[0][0];
    A.y=Floor[0][1];
    A.z=Floor[0][2];
    B.x=Floor[1][0];
    B.y=Floor[1][1];
    B.z=Floor[1][2];
    C.x=Floor[2][0];
    C.y=Floor[2][1];
    C.z=Floor[2][2];
    bresenhamline(Floor[0][0],Floor[0][1],Floor[1][0],Floor[1][1],8);
    bresenhamline(Floor[1][0],Floor[1][1],Floor[2][0],Floor[2][1],8);
    bresenhamline(Floor[2][0],Floor[2][1],Floor[3][0],Floor[3][1],8);
    bresenhamline(Floor[3][0],Floor[3][1],Floor[0][0],Floor[0][1],8);

    for(int i = 0; i<N; i++){
        X.x = light[0];
        X.y = light[1];
        X.z = light[2];
        Y.x = fig[i][0];
        Y.y = fig[i][1];
        Y.z = fig[i][2];
        pt = PlaneIntersectLine(A,B,C,X,Y);
        shadow[i][0]=pt.x;
        shadow[i][1]=pt.y;
        shadow[i][2]=pt.z;
    }
        for(int i = 0; i<N; i++){
        X.x = light[0];
        X.y = light[1];
        X.z = light[2];
        Y.x = fig2[i][0];
        Y.y = fig2[i][1];
        Y.z = fig2[i][2];
        pt = PlaneIntersectLine(A,B,C,X,Y);
        shadow2[i][0]=pt.x;
        shadow2[i][1]=pt.y;
        shadow2[i][2]=pt.z;
    }
    draw_shdow(shadow,shadow2);

}
xyz CreateVector(xyz A, xyz B)
{
xyz m;
m.x = B.x-A.x;
m.y = B.y-A.y;
m.z = B.z-A.z;
return m;
}
xyz VectorProduct(xyz A, xyz B)
{
xyz VP;
VP.x = A.y*B.z-B.y*A.z;
VP.y = A.z*B.x-B.z*A.x;
VP.z = A.x*B.y-B.x*A.y;
return VP;
}
double DotProduct(xyz A, xyz B)
{
double vsp;
vsp = A.x*B.x + A.y*B.y + A.z*B.z;
return vsp;
}
void Normalize(xyz &A)
{
double R,mlr;

mlr = sqrt(pow(A.x,2)+pow(A.y,2)+pow(A.z,2));
A.x = A.x/mlr;
A.y = A.y/mlr;
A.z = A.z/mlr;
}
xyz PlaneIntersectLine(xyz A, xyz B, xyz C, xyz X, xyz Y)
{
xyz rv,Nl,V,W;
double e,d;
NotIntersectPlaneLine = true;
Nl =  VectorProduct(CreateVector(A,B),CreateVector(A,C));
Normalize(Nl);
V = CreateVector(X,A);
d =  DotProduct(Nl,V);
W = CreateVector(X,Y);
e =  DotProduct(Nl,W);
if(e != 0)
{
rv.x = X.x  +  W.x*d/e;
rv.y = X.y  +  W.y*d/e;
rv.z = X.z  +  W.z*d/e;
NotIntersectPlaneLine = false;
}
return rv;
}
void scale(double fig[N][M], double S) {
	double x = 0, y = 0, z = 0;
	x = aver(fig,0);
	y = aver(fig,1);
	z = aver(fig,2);
	double	Sx_Sy[M][M] = { {S,0,0,0},
			  {0,S,0,0},
			  {0,0,S,0},
			  {x*(1 - S),y*(1 - S),z*(1 - S),1} };
	mult(fig, Sx_Sy);
}

void draw(double fig1[N][M],double fig2[N][M]) {
    trianglePaint(fig1[0][0],fig1[0][1],fig1[0][2],fig1[1][0],fig1[1][1],fig1[1][2],fig1[2][0],fig1[2][1],fig1[2][2],5);
	trianglePaint(fig1[1][0],fig1[1][1],fig1[1][2],fig1[4][0],fig1[4][1],fig1[4][2],fig1[2][0],fig1[2][1],fig1[2][2],5);
	trianglePaint(fig1[1][0],fig1[1][1],fig1[1][2],fig1[3][0],fig1[3][1],fig1[3][2],fig1[4][0],fig1[4][1],fig1[4][2],5);
	trianglePaint(fig1[8][0],fig1[8][1],fig1[8][2],fig1[6][0],fig1[6][1],fig1[6][2],fig1[9][0],fig1[9][1],fig1[9][2],5);
	trianglePaint(fig1[6][0],fig1[6][1],fig1[6][2],fig1[7][0],fig1[7][1],fig1[7][2],fig1[9][0],fig1[9][1],fig1[9][2],5);
	trianglePaint(fig1[6][0],fig1[6][1],fig1[6][2],fig1[5][0],fig1[5][1],fig1[5][2],fig1[7][0],fig1[7][1],fig1[7][2],5);
	trianglePaint(fig1[3][0],fig1[3][1],fig1[3][2],fig1[9][0],fig1[9][1],fig1[9][2],fig1[4][0],fig1[4][1],fig1[4][2],11);
	trianglePaint(fig1[3][0],fig1[3][1],fig1[3][2],fig1[8][0],fig1[8][1],fig1[8][2],fig1[9][0],fig1[9][1],fig1[9][2],11);
	trianglePaint(fig1[4][0],fig1[4][1],fig1[4][2],fig1[7][0],fig1[7][1],fig1[7][2],fig1[2][0],fig1[2][1],fig1[2][2],7);
	trianglePaint(fig1[4][0],fig1[4][1],fig1[4][2],fig1[9][0],fig1[9][1],fig1[9][2],fig1[7][0],fig1[7][1],fig1[7][2],7);
	trianglePaint(fig1[2][0],fig1[2][1],fig1[2][2],fig1[5][0],fig1[5][1],fig1[5][2],fig1[0][0],fig1[0][1],fig1[0][2],11);
    trianglePaint(fig1[2][0],fig1[2][1],fig1[2][2],fig1[7][0],fig1[7][1],fig1[7][2],fig1[5][0],fig1[5][1],fig1[5][2],11);
	trianglePaint(fig1[0][0],fig1[0][1],fig1[0][2],fig1[6][0],fig1[6][1],fig1[6][2],fig1[1][0],fig1[1][1],fig1[1][2],7);
	trianglePaint(fig1[0][0],fig1[0][1],fig1[0][2],fig1[5][0],fig1[5][1],fig1[5][2],fig1[6][0],fig1[6][1],fig1[6][2],7);
	trianglePaint(fig1[1][0],fig1[1][1],fig1[1][2],fig1[6][0],fig1[6][1],fig1[6][2],fig1[8][0],fig1[8][1],fig1[8][2],11);
	trianglePaint(fig1[1][0],fig1[1][1],fig1[1][2],fig1[8][0],fig1[8][1],fig1[8][2],fig1[3][0],fig1[3][1],fig1[3][2],11);

	trianglePaint(fig2[0][0],fig2[0][1],fig2[0][2],fig2[1][0],fig2[1][1],fig2[1][2],fig2[2][0],fig2[2][1],fig2[2][2],6);
	trianglePaint(fig2[3][0],fig2[3][1],fig2[3][2],fig2[5][0],fig2[5][1],fig2[5][2],fig2[4][0],fig2[4][1],fig2[4][2],6);
    trianglePaint(fig2[0][0],fig2[0][1],fig2[0][2],fig2[3][0],fig2[3][1],fig2[3][2],fig2[4][0],fig2[4][1],fig2[4][2],10);
    trianglePaint(fig2[0][0],fig2[0][1],fig2[0][2],fig2[4][0],fig2[4][1],fig2[4][2],fig2[1][0],fig2[1][1],fig2[1][2],10);
    trianglePaint(fig2[1][0],fig2[1][1],fig2[1][2],fig2[4][0],fig2[4][1],fig2[4][2],fig2[5][0],fig2[5][1],fig2[5][2],4);
    trianglePaint(fig2[1][0],fig2[1][1],fig2[1][2],fig2[5][0],fig2[5][1],fig2[5][2],fig2[2][0],fig2[2][1],fig2[2][2],4);
    trianglePaint(fig2[2][0],fig2[2][1],fig2[2][2],fig2[5][0],fig2[5][1],fig2[5][2],fig2[3][0],fig2[3][1],fig2[3][2],3);
    trianglePaint(fig2[2][0],fig2[2][1],fig2[2][2],fig2[3][0],fig2[3][1],fig2[3][2],fig2[0][0],fig2[0][1],fig2[0][2],3);
}
void draw_shdow(double fig1[N][M],double fig2[N][M]) {
    trianglePaint(fig1[0][0],fig1[0][1],fig1[0][2],fig1[1][0],fig1[1][1],fig1[1][2],fig1[2][0],fig1[2][1],fig1[2][2],8);
	trianglePaint(fig1[1][0],fig1[1][1],fig1[1][2],fig1[4][0],fig1[4][1],fig1[4][2],fig1[2][0],fig1[2][1],fig1[2][2],8);
	trianglePaint(fig1[1][0],fig1[1][1],fig1[1][2],fig1[3][0],fig1[3][1],fig1[3][2],fig1[4][0],fig1[4][1],fig1[4][2],8);
	trianglePaint(fig1[8][0],fig1[8][1],fig1[8][2],fig1[6][0],fig1[6][1],fig1[6][2],fig1[9][0],fig1[9][1],fig1[9][2],8);
	trianglePaint(fig1[6][0],fig1[6][1],fig1[6][2],fig1[7][0],fig1[7][1],fig1[7][2],fig1[9][0],fig1[9][1],fig1[9][2],8);
	trianglePaint(fig1[6][0],fig1[6][1],fig1[6][2],fig1[5][0],fig1[5][1],fig1[5][2],fig1[7][0],fig1[7][1],fig1[7][2],8);
	trianglePaint(fig1[3][0],fig1[3][1],fig1[3][2],fig1[9][0],fig1[9][1],fig1[9][2],fig1[4][0],fig1[4][1],fig1[4][2],8);
	trianglePaint(fig1[3][0],fig1[3][1],fig1[3][2],fig1[8][0],fig1[8][1],fig1[8][2],fig1[9][0],fig1[9][1],fig1[9][2],8);
	trianglePaint(fig1[4][0],fig1[4][1],fig1[4][2],fig1[7][0],fig1[7][1],fig1[7][2],fig1[2][0],fig1[2][1],fig1[2][2],8);
	trianglePaint(fig1[4][0],fig1[4][1],fig1[4][2],fig1[9][0],fig1[9][1],fig1[9][2],fig1[7][0],fig1[7][1],fig1[7][2],8);
	trianglePaint(fig1[2][0],fig1[2][1],fig1[2][2],fig1[5][0],fig1[5][1],fig1[5][2],fig1[0][0],fig1[0][1],fig1[0][2],8);
    trianglePaint(fig1[2][0],fig1[2][1],fig1[2][2],fig1[7][0],fig1[7][1],fig1[7][2],fig1[5][0],fig1[5][1],fig1[5][2],8);
	trianglePaint(fig1[0][0],fig1[0][1],fig1[0][2],fig1[6][0],fig1[6][1],fig1[6][2],fig1[1][0],fig1[1][1],fig1[1][2],8);
	trianglePaint(fig1[0][0],fig1[0][1],fig1[0][2],fig1[5][0],fig1[5][1],fig1[5][2],fig1[6][0],fig1[6][1],fig1[6][2],8);
	trianglePaint(fig1[1][0],fig1[1][1],fig1[1][2],fig1[6][0],fig1[6][1],fig1[6][2],fig1[8][0],fig1[8][1],fig1[8][2],8);
	trianglePaint(fig1[1][0],fig1[1][1],fig1[1][2],fig1[8][0],fig1[8][1],fig1[8][2],fig1[3][0],fig1[3][1],fig1[3][2],8);

	trianglePaint(fig2[0][0],fig2[0][1],fig2[0][2],fig2[1][0],fig2[1][1],fig2[1][2],fig2[2][0],fig2[2][1],fig2[2][2],8);
	trianglePaint(fig2[3][0],fig2[3][1],fig2[3][2],fig2[5][0],fig2[5][1],fig2[5][2],fig2[4][0],fig2[4][1],fig2[4][2],8);
    trianglePaint(fig2[0][0],fig2[0][1],fig2[0][2],fig2[3][0],fig2[3][1],fig2[3][2],fig2[4][0],fig2[4][1],fig2[4][2],8);
    trianglePaint(fig2[0][0],fig2[0][1],fig2[0][2],fig2[4][0],fig2[4][1],fig2[4][2],fig2[1][0],fig2[1][1],fig2[1][2],8);
    trianglePaint(fig2[1][0],fig2[1][1],fig2[1][2],fig2[4][0],fig2[4][1],fig2[4][2],fig2[5][0],fig2[5][1],fig2[5][2],8);
    trianglePaint(fig2[1][0],fig2[1][1],fig2[1][2],fig2[5][0],fig2[5][1],fig2[5][2],fig2[2][0],fig2[2][1],fig2[2][2],8);
    trianglePaint(fig2[2][0],fig2[2][1],fig2[2][2],fig2[5][0],fig2[5][1],fig2[5][2],fig2[3][0],fig2[3][1],fig2[3][2],8);
    trianglePaint(fig2[2][0],fig2[2][1],fig2[2][2],fig2[3][0],fig2[3][1],fig2[3][2],fig2[0][0],fig2[0][1],fig2[0][2],8);
}
void resetVewport(){
	clearviewport();
	nullZbuf();
}

void bresenhamline(int x0,int y0,int x1,int y1, int color){
    int dx = abs(x1 - x0), sx = x0 < x1 ? 1 : -1;
    int dy = abs(y1 - y0), sy = y0 < y1 ? 1 : -1;
    int err = (dx > dy ? dx : -dy) / 2;
    int e2 = err;
    for (;;) {
        putpixel(x0, y0, color);
        if (x0 == x1 && y0 == y1) break;
        e2 = err;
        if (e2 > -dx){
            err -= dy;
            x0 += sx;
        }
        if (e2 < dy) {
            err += dx;
            y0 += sy;
        }
    }
}
void rotateX(double fig[N][M], double angle) {
	double y = 0, z = 0;
	y = aver(fig,1);
	z = aver(fig,2);
	double Angle[M][M] = { {1,0, 0, 0},
			{0 , cos(angle), sin(angle),0},
			{0, -sin(angle), cos(angle), 0},
			{0, y*(1 - cos(angle)) + z * sin(angle), z*(1 - cos(angle)) - y * sin(angle), 1} };
	mult(fig, Angle);
}
void rotateY(double fig[N][M], double angle) {
	double x = 0, y = 0, z = 0;
	x = aver(fig,0);
	z = aver(fig,2);
	double Angle[M][M] = { {cos(angle),0, -sin(angle), 0},
			{0, 1, 0,0},
			{sin(angle), 0, cos(angle), 0},
			{x*(1 - cos(angle)) - z * sin(angle), 0, z*(1 - cos(angle)) + x * sin(angle), 1} };

	mult(fig, Angle);
}

void rotateZ(double fig[N][M], double angle) {
	double x = 0, y = 0;
	x=aver(fig, 0);
	y=aver(fig, 1);
	double Angle[M][M] = { {cos(angle), sin(angle), 0, 0},
		{ -sin(angle), cos(angle), 0, 0},
			 {0, 0, 1, 0},
			 {x*(1 - cos(angle)) + y * sin(angle), y*(1 - cos(angle)) - x * sin(angle), 0, 1} };
	mult(fig, Angle);
}

double aver(double fig[N][M], int cnt){
	double average=0;
	for (int i = 0; i< N; i++) {
		average += fig[i][cnt];
	}
	return average/N;
}
void nullZbuf(){
    int i=0,j=0;
    for(i=0;i<ScreenWidth;i++)
        for(j=0;j<ScreenHeigh;j++)
            Zbuffer[i][j]=1000.0;
}
void trianglePaint(int x1,int y1,int z1,int x2,int y2,int z2, int x3,int y3,int z3, int color){
    int Ymax;
    int Ymin;
    int Xmax;
    int Xmin;
    int a,b,c,d,z;
    Ymax = y1<y2 ? y2:y1;
    Ymax = Ymax<y3 ? y3:Ymax;
    Ymin = y1<y2 ? y1:y2;
    Ymin = Ymin<y3 ? Ymin:y3;
    Xmax = x1<x2 ? x2:x1;
    Xmax = Xmax<x3 ? x3:Xmax;
    Xmin = x1<x2 ? x1:x2;
    Xmin = Xmin<x3 ? Xmin:x3;
    Vertex v0 = {x1,y1};
    Vertex v1 = {x2,y2};
    Vertex v2 = {x3,y3};

    a = ((y2-y1)*(z3-z1)-(y3-y1)*(z2-z1));
    b = ((x3-x1)*(z2-z1)-(x2-x1)*(z3-z1));
    c = ((x2-x1)*(y3-y1)-(y2-y1)*(x3-x1));
    d = -x1*((y2-y1)*(z3-z1)-(y3-y1)*(z2-z1))-y1*((x3-x1)*(z2-z1)-(x2-x1)*(z3-z1))-z1*((x2-x1)*(y3-y1)-(y2-y1)*(x3-x1));
    for (unsigned int y = Ymin; y < Ymax; ++y)
    {
        for (unsigned int x = Xmin; x < Xmax; ++x)
        {
            Point p = {(float)x, (float)y};
            float e10 = Edge(v1.position, v0.position, p);
            float e21 = Edge(v2.position, v1.position, p);
            float e02 = Edge(v0.position, v2.position, p);
            if (e10 >= 0.0f && e21 >= 0.0f && e02 >= 0.0f)
            {
                if (c!=0){
                    z = (-(a*x+b*y+d))/c;
                    if(Zbuffer[x][y]>z)
                    {
                        putpixel(x, y, color);
                        Zbuffer[x][y] = z;
                    }
                }
            }
        }
    }
}
