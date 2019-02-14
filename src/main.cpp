
const int bottomBorder = 600;
const int rightBorder = 600;
const int numofparticles = 100;
const int h = 30;
const double pi = 3.14159265359f;
const double e = 2.71828182846f;
const double g = -9.80665f;
const double ny = 0;//.000001827f;
const double mass = 1000.0f;
const double dens = 3000.0f;
const double const_dv = 0.5;
const double const_press = 0.5;

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <SFML\Graphics.hpp>
#include <vector>
#include <omp.h>
#include <amp.h>

using namespace std;
using namespace concurrency;


struct Vector3 {

	double x, y, z;

	Vector3(double xx, double yy, double zz)
	{
		x = xx;
		y = yy;
		z = zz;
	};
	Vector3(const Vector3& vec)
	{
		x = vec.x;
		y = vec.y;
		z = vec.z;
	};
	Vector3()
	{
		x = 0.0f;
		y = 0.0f;
		z = 0.0f;
	};

	Vector3 operator +  (const Vector3& vec)
	{
		Vector3 result{ 0.0, 0.0, 0.0 };

		result.x = x + vec.x;
		result.y = y + vec.y;
		result.z = z + vec.z;

		return result;
	};

	Vector3 operator * (const double& con)
	{
		Vector3 result;

		result.x = x*con;
		result.y = y*con;
		result.z = z*con;

		return result;
	};

	Vector3 operator / (const double& con)
	{
		Vector3 result;

		result.x = x / con;
		result.y = y / con;
		result.z = z / con;

		return result;
	};

	Vector3 operator - (const Vector3& vec)
	{
		Vector3 result{ 0.0, 0.0, 0.0 };

		result.x = x - vec.x;
		result.y = y - vec.y;
		result.z = z - vec.z;

		return result;
	};

	void print()
	{
		cout << "x= " << x << endl;
		cout << "y= " << y << endl;
		cout << "z= " << z << endl;
		cout << "\n" << endl;
	};

	double length_one()
	{
		return sqrt(powf(x, 2) + powf(y, 2) + powf(z, 2));
	};

	Vector3& normal()
	{
		this->x = x / sqrt(x*x + y*y + z*z);
		this->y = y / sqrt(x*x + y*y + z*z);
		this->z = z / sqrt(x*x + y*y + z*z);

		return *this;
	};
};

const Vector3 G{ 0, 0, -9.80665f };

struct particle {
	Vector3 position;
	Vector3 velocity;
	double pressure;
	double density;
	particle()
	{
		position = { 0.0f, 0.0f, 0.0f };
		velocity = { 0.0f, 0.0f, 0.0f };
		pressure = 10;
		density = 100;
	};

	particle(Vector3 posit, Vector3 vel, double press, double dens, double mas)
	{
		position = posit;
		velocity = vel;
		pressure = press;
		density = dens;
	};
};

double scalar(Vector3 vec1, Vector3 vec2) 
{
	return vec1.x*vec2.x + vec1.y*vec2.y + vec1.z*vec2.z;
}

double length(Vector3 vec1, Vector3 vec2) 
{
	return sqrt(powf(abs(vec1.x - vec2.x), 2) + powf(abs(vec1.y - vec2.y), 2) + powf(abs(vec1.z - vec2.z), 2));
}


Vector3& W_press(Vector3& r)            
{
	Vector3 vec{ 0.0f, 0.0f, 0.0f };

	if (r.length_one() <= h && r.length_one() != 0)
	{
		vec = r.normal()*((45 / (pi*powf(h, 6)))*powf(h - r.length_one(), 3));

		return vec;
	}

	else return vec;
}

double W_vis(Vector3& r)                  
{
	if (r.length_one() <= h )
	{
		return (45 / (pi*pow(h, 6)))*(h - r.length_one());
	}
	else return 0;
}

double W(Vector3& r)												
{																
	if (r.length_one() <= h && r.length_one() != 0)
	{
		return (315 / (64 * pi*powf(h, 9)))*powf(h*h - pow(r.length_one(), 2), 3);
	}
	else return 0;
}

double density(vector<particle>& vect_part, particle& part)           
{
	double sum = 1000;
	for (int i = 0; i < vect_part.size(); ++i)
	{
		sum += W(part.position - vect_part[i].position);
	}
	return sum;
}

double pressure(particle& part)
{
	return 100000 + abs(const_press*(part.position.z)*(part.density - 1000));
}

Vector3 Navier_Stokes_i(vector<particle>& vect_part, particle& part) 
{
	Vector3 vec={ 0.0f, 0.0f, 0.0f }, vec1 = { 0.0f, 0.0f, 0.0f }, vec2={ 0.0f, 0.0f, 0.0f };
#pragma omp parallel
	{
		#pragma omp for
		for (int i = 0; i < vect_part.size(); ++i)
		{
			vect_part[i].density = density(vect_part, vect_part[i]);
			vect_part[i].pressure = pressure(vect_part[i]);
		}
		#pragma omp for
		for (int i = 0; i < vect_part.size(); ++i)
		{
			vec2 = vec2 + W_press(part.position - vect_part[i].position) * ((vect_part[i].pressure + part.pressure) / (2 * vect_part[i].density));
			vec1 = vec1 + (vect_part[i].velocity - part.velocity) * (ny * W(vect_part[i].position - part.position) / (2 * vect_part[i].density));
		}
	};
	vec = (vec2)*(-1) + vec1;
	return vec;
}



int main()
{
	vector<particle> mol;
	vector<Vector3> dv(numofparticles);
	particle a;
	for (int i = 0; i < numofparticles; i++)
	{
		a.position.x = (double)rightBorder * i / RAND_MAX;
		a.position.y = 0;
		a.position.z = -(double)bottomBorder * i / RAND_MAX;
		a.velocity.x = 1.5 - 1.0 * rand() / RAND_MAX;
		a.velocity.z = -0.5 - 1.0 * rand() / RAND_MAX;
		mol.push_back(a);
	};
	sf::CircleShape circle(7);
	circle.setFillColor(sf::Color(50, 50, 200,155));
	sf::RenderWindow window(sf::VideoMode(rightBorder, bottomBorder), "Liquids-SPH v2.0.2");
	double time = 0;
	double dt = 0;
	while (window.isOpen())
	{
		if (sf::Keyboard::isKeyPressed(sf::Keyboard::P))
		{
			system("pause");
		}
		sf::Event event;
		while (window.pollEvent(event))
		{
			if (event.type == sf::Event::Closed)
				window.close();
		}
		dt = 0.5;
		window.clear(sf::Color::White);
		#pragma omp parallel
		{
			#pragma omp for
			for (int i = 0; i < numofparticles; i++)
			{
				sf::Vector2f pos;
				pos.x = (float)mol[i].position.x;
				pos.y = -(float)mol[i].position.z;
				circle.setPosition(pos);
				window.draw(circle);
				if (mol[i].position.x + mol[i].velocity.x * dt * 100 > rightBorder -7|| mol[i].position.x + mol[i].velocity.x * dt * 100 < 7) mol[i].velocity.x = mol[i].velocity.x * (-0.3);
				if (mol[i].position.z + mol[i].velocity.z * dt * 100 < -bottomBorder +7 || mol[i].position.z + mol[i].velocity.z * dt * 100 > -7) mol[i].velocity.z = mol[i].velocity.z * (-0.3);
				#pragma omp for
				for (int k = i; k < numofparticles; ++k)
				{
					if (length(mol[k].position + mol[k].velocity * dt * 100, mol[i].position + mol[i].velocity* dt * 100) <= 5 && length(mol[k].position, mol[i].position) > 5)
					{
						Vector3 vec = { 0.0f,0.0f,0.0f };
						vec = (mol[i].velocity + mol[k].velocity)*0.5;
						mol[i].velocity = vec;
						mol[k].velocity = vec;
					}
				}
				mol[i].position = mol[i].position + mol[i].velocity * dt * 100;
				dv[i] = Navier_Stokes_i(mol, mol[i]) * dt *const_dv;		
			}
			#pragma omp for
			for (int i = 0; i < numofparticles; i++)
			{
				mol[i].velocity = mol[i].velocity + dv[i];
			}	
		}
		window.display();
	}
	system("pause");
	return 0;
}

