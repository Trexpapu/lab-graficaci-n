// rt: un lanzador de rayos minimalista
 // g++ -O3 -fopenmp rt.cpp -o rt
#include <math.h>
#include <stdlib.h>
#include <stdio.h>  
#include <omp.h>
#include <cmath>  // Para la función pow
#include <unistd.h> // Para la función sleep
#include <random>
#include <iostream>


float GETNEXTRAND() {
    // Crea un generador de números aleatorios
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> dist(0.0, 1.0);  // Distribución entre 0 y 1

    return dist(gen);
}

class Vector 
{
public:        
	double x, y, z; // coordenadas x,y,z 
  
	// Constructor del vector, parametros por default en cero
	Vector(double x_= 0, double y_= 0, double z_= 0){ x=x_; y=y_; z=z_; }
  
	// operador para suma y resta de vectores
	Vector operator+(const Vector &b) const { return Vector(x + b.x, y + b.y, z + b.z); }
	Vector operator-(const Vector &b) const { return Vector(x - b.x, y - b.y, z - b.z); }
	// operator multiplicacion vector y escalar 
	Vector operator*(double b) const { return Vector(x * b, y * b, z * b); }
  
	// operator % para producto cruz
	Vector operator%(const Vector&b)const {return Vector(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x);}
	
	// producto punto con vector b
	double dot(const Vector &b) const { return x * b.x + y * b.y + z * b.z; }

	// producto elemento a elemento (Hadamard product)
	Vector mult(const Vector &b) const { return Vector(x * b.x, y * b.y, z * b.z); }
	
	// normalizar vector 
	Vector& normalize(){ return *this = *this * (1.0 / sqrt(x * x + y * y + z * z)); }
};
typedef Vector Point;
typedef Vector Color;

class Ray 
{ 
public:
	Point o;
	Vector d; // origen y direcccion del rayo
	Ray(Point o_, Vector d_) : o(o_), d(d_) {} // constructor
};

double menor(double a, double b){
	if(a < b){
		return a;
	}else{
		return b;
	}

	return 0;
}

void coordinateSystem(const Vector &n, Vector &s, Vector &t) {
	if (std::abs(n.x) > std::abs(n.y)) {
		float invLen = 1.0f / std::sqrt(n.x * n.x + n.z * n.z);
		t = Vector(n.z * invLen, 0.0f, -n.x * invLen);
	}

	else 
	{
		float invLen = 1.0f / std::sqrt(n.y * n.y + n.z * n.z);
		t = Vector(0.0f, n.z * invLen, -n.y * invLen);		
	}
	s = n % t;
}

Vector LocalesAGlobales(const Vector &n,Vector &s,Vector &t,const Vector &local) {
	Vector global;
	global.x = s.x * local.x + t.x * local.y + n.x * local.z;
	global.y = s.y * local.x + t.y * local.y + n.y * local.z;
	global.z = s.z * local.x + t.z * local.y + n.z * local.z;
	return global;
}

Vector GlobalesALocales(const Vector &n,Vector &s,Vector &t,const Vector &global) {
	Vector local;
	local.x = s.dot(global);
	local.y = t.dot(global);
	local.z = n.dot(global);
	return local;
}


class Sphere 
{
public:
	double r;	// radio de la esfera
	Point p;	// posicion
	Color c;	// color  
	Color emluz; // color de la luz

	Sphere(double r_, Point p_, Color c_, Color emluz_): r(r_), p(p_), c(c_), emluz(emluz_) {} 

  
	// PROYECTO 1
	// determina si el rayo intersecta a esta esfera
	double intersect(const Ray &ray) const {
		// regresar distancia si hay intersección
		// regresar 0.0 si no hay interseccion

		//primero calculamos el discriminante para saber si vale la pena o no seguir

		double discriminante = sqrt(
			pow((ray.o - p).dot(ray.d), 2) // primera parte de la raiz cuadrada
			-
			(ray.o  - p).dot(ray.o - p) //segunda parte de la raiz cuadrada
			+ 
			pow(r, 2) //tercera parte de la raiz cuadrada
		);

		if(discriminante >= 0){//entonces si hay un choque  o intersecciòn con la esfera
			//calculamos la parte que no es discriminante
			double no_discriminante = - ((ray.o - p).dot(ray.d));

			//primero el valor con la operacion +
			double t1 = no_discriminante + discriminante;
			//segundo el valor con la operacion -
			double t2 = no_discriminante - discriminante;
            

            
            if(t1 > 0.0 && t2 > 0.0){//en el caso de que ambas sean t sean positivas retornamos el mas pequeño
                return menor(t1, t2);
            }else if(t1 > 0.0){
                return t1;
            }else if(t2 > 0.0){
                return t2;
            }

		}

		return 0.0;
	}
};

Sphere spheres[] = {
	//Escena: radio, posicion, color 
	Sphere(1e5,  Point(-1e5 - 49, 0, 0),   Color(.75, .25, .25), Color()), // pared izq
	Sphere(1e5,  Point(1e5 + 49, 0, 0),    Color(.25, .25, .75), Color()), // pared der
	Sphere(1e5,  Point(0, 0, -1e5 - 81.6), Color(.25, .75, .25), Color()), // pared detras
	Sphere(1e5,  Point(0, -1e5 - 40.8, 0), Color(.25, .75, .75), Color()), // suelo
	Sphere(1e5,  Point(0, 1e5 + 40.8, 0),  Color(.75, .75, .25), Color()), // techo
	Sphere(16.5, Point(-23, -24.3, -34.6), Color(.2, .3, .4), Color()), // esfera abajo-izq
	Sphere(16.5, Point(23, -24.3, -3.6),   Color(.4, .3, .2), Color()), // esfera abajo-der
	Sphere(10.5, Point(0, 24.3, 0),        Color(1, 1, 1), Color(10, 10, 10)) // arriba
};

// limita el valor de x a [0,1]
inline double clamp(const double x) { 
	if(x < 0.0)
		return 0.0;
	else if(x > 1.0)
		return 1.0;
	return x;
}

// convierte un valor de color en [0,1] a un entero en [0,255]
inline int toDisplayValue(const double x) {
	return int( pow( clamp(x), 1.0/2.2 ) * 255 + .5); 
}

// PROYECTO 1
// calcular la intersección del rayo r con todas las esferas
// regresar true si hubo una intersección, falso de otro modo
// almacenar en t la distancia sobre el rayo en que sucede la interseccion
// almacenar en id el indice de spheres[] de la esfera cuya interseccion es mas cercana
inline bool intersect(const Ray &r, double &t, int &id) {
	double inf = 1e20; //valor grande que representa el infinito
	t = inf;
	bool hit = false;
	for(int i = 0; i < sizeof(spheres) / sizeof(Sphere); i++){
		double distance = spheres[i].intersect(r);
		if ( distance > 1e-8 && distance < t){
			t = distance; //guarda la distancia
			id = i; //guarda el indice de la esfera
			hit = true;
		}
	}
	return hit;
}




class MonteCarlo{
	Ray ray;
	Point x; 
	Vector n;
	Sphere obj;

	public:

	MonteCarlo(Ray ray_, Point x_, Vector n_): ray(ray_), x(x_), n(n_), obj(0, Point(), Color(), Color()){}

	Vector uniformeEsferico(Sphere sphere) {
		
		double thetaj = acos(1 - 2 * GETNEXTRAND()); 
		double pihemisj = 2 * M_PI * GETNEXTRAND(); 
		double pwj = 1 / (4 * M_PI); 

		Vector wi (
			cos(pihemisj) * sin(thetaj),
			sin(pihemisj) * sin(thetaj),
			cos(thetaj)
		);

		Vector s, t;
		coordinateSystem(n, s, t);

		Vector global = LocalesAGlobales(n, s, t, wi);

		double costhetai = n.dot(global);
		Color fr = sphere.c * (1 / M_PI);
		Color cle = le(global);
		Vector result = ((fr.mult(cle) * costhetai)) * (1 / pwj);
		return result;
	}

	Vector uniformeHemisferico(Sphere sphere) {
		
		double thetaj = acos(GETNEXTRAND()); 
		double pihemisj = 2 * M_PI * GETNEXTRAND(); 
		double pwj = 1 / (2 * M_PI); 

		Vector wi (
			cos(pihemisj) * sin(thetaj),
			sin(pihemisj) * sin(thetaj),
			cos(thetaj)
		);

		Vector s, t;
		coordinateSystem(n, s, t);

		Vector global = LocalesAGlobales(n, s, t, wi);

		double costhetai = n.dot(global);
		Color fr = sphere.c * (1 / M_PI);
		Color cle = le(global);
		Vector result = ((fr.mult(cle) * costhetai)) * (1 / pwj);
		return result;
	}

	Vector cosenoHemisferico(Sphere sphere) {
		double thetaj = acos(sqrt(1 - GETNEXTRAND()));
		double pihemisj = 2 * M_PI * GETNEXTRAND();
		double pwj = ((1/(M_PI * cos(thetaj)))); 
		
		Vector wi (
			cos(pihemisj) * sin(thetaj),
			sin(pihemisj) * sin(thetaj),
			cos(thetaj)
		);

		Vector s, t;
		coordinateSystem(n, s, t);

		Vector global = LocalesAGlobales(n, s, t, wi);

		double costhetai = n.dot(global);
		Color fr = sphere.c * (1 / M_PI);
		Color cle = le(global);
		Vector result = ((fr.mult(cle) * costhetai)) * (1 / pwj);
		return result;
	}


	private: 
	Color le(const Vector global){
		double t;
		int id = 0;
		Color cle = Color();
		Ray rayo = Ray(x, global);
		if(intersect(rayo, t, id)){
			cle = spheres[id].emluz;
		}
		return cle;
	}
};

// Calcula el valor de color para el rayo dado
Color shade(const Ray &r) {
	double t;
	int id = 0;
	// determinar que esfera (id) y a que distancia (t) el rayo intersecta
	if (!intersect(r, t, id))
		return Color();	// el rayo no intersecto objeto, return Vector() == negro
  
	const Sphere &obj = spheres[id];
	

	if (obj.emluz.x > 0 || obj.emluz.y > 0 || obj.emluz.z > 0) {
		return obj.emluz;
	}

	// PROYECTO 1
	
	Point x = r.o + r.d * t;

	
	Vector n = (x - obj.p).normalize();

	MonteCarlo monteCarlo(r, x, n);

	Color colorValue = monteCarlo.uniformeEsferico(obj) + obj.emluz;

	//Color colorValue = monteCarlo.uniformeHemisferico(obj) + obj.emluz;

	//Color colorValue = monteCarlo.cosenoHemisferico(obj) + obj.emluz;

	return colorValue; 
}


int main(int argc, char *argv[]) {
	

	int muestra = 32; // número de muestras

	int w = 1024, h = 768; // image resolution
  
	// fija la posicion de la camara y la dirección en que mira
	Ray camera( Point(0, 11.2, 214), Vector(0, -0.042612, -1).normalize() );

	// parametros de la camara
	Vector cx = Vector( w * 0.5095 / h, 0., 0.); 
	Vector cy = (cx % camera.d).normalize() * 0.5095;
  
	// auxiliar para valor de pixel y matriz para almacenar la imagen
	Color *pixelColors = new Color[w * h];

	// PROYECTO 1
	// usar openmp para paralelizar el ciclo: cada hilo computara un renglon (ciclo interior),
	int Porcentaje = 0; // variable para ver cuántas filas ya fueron procesadas
	#pragma omp parallel for schedule(dynamic, 1)
	for(int y = 0; y < h; y++) 
	{ 
		// recorre todos los pixeles de la imagen
		// Actualizar el progreso 
		#pragma omp atomic
		Porcentaje++;

		
		#pragma omp critical
		{
			fprintf(stderr, "\r%5.2f%%", 100.0 * Porcentaje / h);
		}
		
		for(int x = 0; x < w; x++ ) {
			int idx = (h - y - 1) * w + x; 
			Color pixelValue = Color(); 
			for (int mt = 0; mt < muestra; mt++) {
				Vector cameraRayDir = cx * (double(x) / w - 0.5) + cy * (double(y) / h - 0.5) + camera.d;

				// llamar a shade para obtener el color del rayo
				pixelValue = pixelValue + shade(Ray(camera.o, cameraRayDir.normalize()));
			}
			pixelValue = pixelValue * (1.0 / muestra);

			// limitar los tres valores de color del pixel a [0,1]
			pixelColors[idx] = Color(clamp(pixelValue.x), clamp(pixelValue.y), clamp(pixelValue.z));
		}
	}

	fprintf(stderr,"\n");

	// PROYECTO 1
	// Investigar formato ppm
	FILE *f = fopen("image.ppm", "w");
	// escribe cabecera del archivo ppm, ancho, alto y valor maximo de color
	fprintf(f, "P3\n%d %d\n%d\n", w, h, 255); 
	for (int p = 0; p < w * h; p++) 
	{ // escribe todos los valores de los pixeles
			fprintf(f,"%d %d %d ", toDisplayValue(pixelColors[p].x), toDisplayValue(pixelColors[p].y), 
				toDisplayValue(pixelColors[p].z));
	}
	fclose(f);

	delete[] pixelColors;

	

	return 0;
}