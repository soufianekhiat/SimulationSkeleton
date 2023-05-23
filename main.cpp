#include <SFML/Graphics.hpp>
#include <vector>
#include <thread>
#include <random>
#include <cmath>
#include <iostream>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <imgui.h>
#include <imgui-SFML.h>

#define M_PIf ((float)M_PI)

constexpr float dt = 0.0125f;

constexpr int K_PARTICLES = 300;
constexpr int WORLD_WIDTH = 1920;
constexpr int WORLD_HEIGTH = 1080;
constexpr int DOT_SIZE = 2;
constexpr float DOT_SIZEf = static_cast< float >( DOT_SIZE );
static sf::Vector2f DOT_OFSET = sf::Vector2f(DOT_SIZEf, DOT_SIZEf);

float sigmaAttraction = 5.0f;

// Simple struct to hold particle data
struct Particle {
	sf::Vector2f position;
	sf::Vector2f previousPosition;
	sf::Vector2f velocity; // For initialization only
	float orientation;
	float previousOrientatation;
	//float dorientation_dt;
	float ddorientation_dt2;
	//float angularVelocity;
	sf::Vector2f force;
	//float torque;
	int type;
	//std::vector<Particle*> linked;
	//Having view relating object in the model object is not ideal
	//but in SFML not rebuilding the graphical object is significantly faster
	sf::CircleShape shape;
};

__forceinline
sf::Vector2f multiply(sf::Vector2f const& vec, float const scalar) {
	return sf::Vector2f(vec.x * scalar, vec.y * scalar);
}

__forceinline
float angleFromVector(sf::Vector2f v) {
	return std::atan2(v.y, v.x);
}

__forceinline
sf::Vector2f unitVectorFromAngle(float iAngle) {
	return sf::Vector2f(std::cos(iAngle), std::sin(iAngle));
}

__forceinline
float dot(sf::Vector2f const& a, sf::Vector2f const& b) {
	// == cos(angle(a, b))
	return a.x * b.x + a.y * b.y;
}

__forceinline
float normSqr(sf::Vector2f const& vec) {
	return dot(vec, vec);
}

__forceinline
float norm(sf::Vector2f const& vec) {
	return std::sqrt(normSqr(vec));
}

void normalize(sf::Vector2f& vec)
{
	float const normVec = normSqr(vec);
	if ( std::abs(normVec) > 1e-8f )
		vec /= std::sqrt(normVec);
}

sf::Vector2f normalized(sf::Vector2f const& _vec)
{
	sf::Vector2f v = _vec;
	float const normVec = normSqr(v);
	if ( std::abs(normVec) > 1e-8f )
		v /= std::sqrt(normVec);

	return v;
}

__forceinline
float saturate(float x)
{
	return x < 0.0f ? 0.0f : ( x > 1.0f ? 1.0f : x );
}

sf::Color getColor(int type) {
	// Convert the type to a hue value between 0 and 360 degrees
	float hue = ( ( float )( type % 16 ) ) * ( 360.0f / 16.0f );

	// For simplicity, we'll keep saturation and lightness constant
	float saturation = 1.0f;
	float lightness = 0.5f;

	// Convert from HSL to RGB using the formula
	float c = ( 1.0f - std::abs(2.0f * lightness - 1.0f) ) * saturation;
	float x = c * ( 1.0f - std::abs(std::fmod(hue / 60.0f, 2.0f) - 1.0f) );
	float m = lightness - c / 2.0f;

	float r = 0.0f, g = 0.0f, b = 0.0f;
	if ( 0.0f <= hue && hue < 60.0f ) {
		r = c, g = x, b = 0.0f;
	}
	else if ( 60.0f <= hue && hue < 120.0f ) {
		r = x, g = c, b = 0.0f;
	}
	else if ( 120.0f <= hue && hue < 180.0f ) {
		r = 0.0f, g = c, b = x;
	}
	else if ( 180.0f <= hue && hue < 240.0f ) {
		r = 0.0f, g = x, b = c;
	}
	else if ( 240.0f <= hue && hue < 300.0f ) {
		r = x, g = 0.0f, b = c;
	}
	else if ( 300.0f <= hue && hue < 360.0f ) {
		r = c, g = 0.0f, b = x;
	}

	return sf::Color(( sf::Uint8 )std::round(saturate(r + m) * 255),
					 ( sf::Uint8 )std::round(saturate(g + m) * 255),
					 ( sf::Uint8 )std::round(saturate(b + m) * 255));
}

float KernelCubicSpline(float const distance, float const smoothDistance)
{
	float const q = distance / smoothDistance;

	if ( 0.0f <= q || q <= 0.5f )
	{
		return 6.0f * ( q * ( q * q - q ) ) + 1.0f;
	}
	else if ( 0.5f < q || q <= 1.0f )
	{
		float const _q = 1.0f - q;
		return 2.0f * _q * _q * _q;
	}
	else
	{
		return 0.0f;
	}
}

struct Model {
	std::vector<Particle> particles;
	std::random_device rd;
	std::mt19937_64 gen;

	Model()
	:	gen(rd())
	{
		init();
	}

	void init()
	{
		// Initialize random number generator
		std::uniform_int_distribution<int> disType(0, 15);
		std::uniform_real_distribution<float> disX(-( ( float )WORLD_WIDTH ) / 16.0f, ( ( float )WORLD_WIDTH ) / 16.0f);
		std::uniform_real_distribution<float> disY(-( ( float )WORLD_HEIGTH ) / 16.0f, ( ( float )WORLD_HEIGTH ) / 16.0f);

		std::uniform_real_distribution<float> dis01(-0.0f, 1.0f);
		std::uniform_real_distribution<float> dis_11(-1.0f, 1.0f);

		// Initialize particles
		particles.clear();
		for ( int i = 0; i < K_PARTICLES; ++i ) {
			Particle p;
			p.position = sf::Vector2f(disX(gen), disY(gen));
			p.velocity = normalized(sf::Vector2f(2.0f * dis_11(gen), 2.0f * dis_11(gen))) * 0.1f;
			p.previousPosition = p.position - p.velocity * dt;
			//p.velocity = normalized(sf::Vector2f(2.0f * dis_11(gen), 2.0f * dis_11(gen))) * 0.1f;
			p.orientation = dis01(gen) * 2.0f * M_PIf;
			p.previousOrientatation = p.orientation;
			//p.dorientation_dt = dis_11(gen) * 0.01f;
			p.ddorientation_dt2 = 0.0f;
			//p.angularVelocity = 0.0;
			p.type = disType(gen);
			p.shape = sf::CircleShape(DOT_SIZEf);
			p.shape.setOrigin(DOT_SIZEf, DOT_SIZEf);
			p.shape.setPosition(p.position);
			particles.push_back(p);
		}
	}

	void step()
	{
		//const float attractionStrength = 0.05f;
		//constexpr float interactionRadius = 30.0f;
		constexpr float interactionRadius = (2.0f * DOT_SIZEf) * 5.0f;
		constexpr float interactionRadiusSqr = interactionRadius * interactionRadius;
		//const float attractionMinThreshold = 5.0f;
		//const float linkingThreshold = 30.0f;

		constexpr float rSize = 2.0f * DOT_SIZEf;
		constexpr float rSizeSqr = rSize * rSize;

		// Calculate the force and torque on particle p due to all other particles
		for ( auto& p : particles )
		{
			p.force = sf::Vector2f(0.0f, 0.0f);
			p.ddorientation_dt2 = 0.0f;
			for ( Particle& other : particles )
			{
				if ( &other != &p )
				{
					// Avoid self-interaction
					sf::Vector2f r = other.position - p.position;
					float const rNormSqr = normSqr(r);
					// Consider only particles within the interaction radius
					if ( rNormSqr < interactionRadiusSqr )
					{
						// Calculate force and torque using appropriate model
						//sf::Vector2f force(0.0f, 0.0f);
						//float torque = 0.0f;
						//p.force += force;

						float const dist = std::sqrt(rNormSqr);
						//p.force += intensity * r * ( -100.0f / dist );
						r /= dist;

						// Repulsive force
						if ( rNormSqr < rSizeSqr )
						{
							sf::Vector2f newPos = p.position - 0.125f * DOT_SIZEf * r;
							p.previousPosition = p.position;
							p.position = newPos;
							p.force = sf::Vector2f(p.force.x * 0.125f, p.force.y * 0.125f);
						}
						else
						{
							float const intensity = KernelCubicSpline(dist, interactionRadius);
							p.force += 0.125f * intensity * r;
						}

						p.ddorientation_dt2 -= 0.0125f * saturate(dot(unitVectorFromAngle(p.orientation), unitVectorFromAngle(other.orientation)));
					}
				}
			}

		//	// Containing forces
		//	if ( p.position.x > WORLD_WIDTH / 2 ) p.force += sf::Vector2f(-0.1f, 0.0f);
		//	if ( p.position.x < -WORLD_WIDTH / 2 ) p.force += sf::Vector2f(0.1f, 0.0f);
		//	if ( p.position.y > WORLD_HEIGTH / 2 ) p.force += sf::Vector2f(0.0f, -0.1f);
		//	if ( p.position.y < -WORLD_HEIGTH / 2 ) p.force += sf::Vector2f(0.0f, 0.1f);
		}

		for ( auto& p : particles ) {
			//sf::Vector2f acceleration = p.force; // / mass (mass == 1)

			// Symplectic Euler 1:
			//p.velocity += acceleration * dt;
			//p.position += p.velocity * dt;
			//
			//p.dorientation_dt += p.ddorientation_dt2 * dt;
			//p.orientation += p.dorientation_dt * dt;

			// Verlet:
			sf::Vector2f const position = 2.0f * p.position - p.previousPosition + p.force * dt * dt;
			float const orientation = 2.0f * p.orientation - p.previousOrientatation + p.ddorientation_dt2 * dt * dt;

			p.previousPosition = p.position;
			p.position = position;
			p.previousOrientatation = p.orientation;
			p.orientation = orientation;

			p.shape.setPosition(p.position);
			p.shape.setFillColor(getColor(p.type));
		}

		//for ( auto& p : particles ) {
		//	// Calculate the acceleration and angular acceleration (assuming mass and moment of inertia = 1)
		//	sf::Vector2f acceleration = p.force;
		//	//float angularAcceleration = p.torque;
		//	float dt = 0.01f;
		//
		//	// Using naive algo
		//	p.velocity = p.velocity + acceleration * dt;
		//	float maxVelocity = 50.0f;
		//	if ( norm(p.velocity) > maxVelocity ) {
		//		normalize(p.velocity);
		//		p.velocity *= maxVelocity;
		//	}
		//	//Force attract to center to incentive interactions
		//	p.velocity -= multiply(p.position, 0.01f);
		//	//Sticky dissipative space and other limits
		//	p.velocity = multiply(p.velocity, 0.90f);
		//	float maxAngularVelocity = 10.0f;
		//	p.position = p.position + p.velocity * dt;
		//	//p.orientation = p.orientation + p.angularVelocity * dt;
		//
		//	// Update the particle's shape position
		//	p.shape.setPosition(p.position);
		//	p.shape.setFillColor(getColor(p.type));
		//}
	}
};

void drawModel(sf::RenderWindow& ioWindow, Model const& iModel) {
	for ( auto const& p : iModel.particles )
	{
		/*for (auto& other : p.linked) {
			sf::VertexArray lines(sf::Lines, 2);
			lines[0].position = p.position + DOT_OFSET;
			lines[1].position = other->position + DOT_OFSET;
			lines[0].color = sf::Color::White;
			lines[1].color = sf::Color::White;
			ioWindow.draw(lines);
		}*/

		float tailLength = 10.0f;
		sf::Vertex line[] =
		{
			sf::Vertex(p.position, sf::Color::White),
			sf::Vertex(p.position + tailLength * unitVectorFromAngle(p.orientation), sf::Color::White)
		};

		//Force debug
		//sf::Vertex lineF[] =
		//{
		//    sf::Vertex(p.position, sf::Color::Red),
		//    sf::Vertex(p.position + p.force, sf::Color::Red)
		//};
		//ioWindow.draw(lineF, 2, sf::Lines);

		//Speed debug
		//sf::Vertex lineV[] =
		//{
		//    sf::Vertex(p.position, sf::Color::Green),
		//    sf::Vertex(p.position + 10.0f*p.velocity, sf::Color::Green)
		//};
		//ioWindow.draw(lineV, 2, sf::Lines);

		ioWindow.draw(line, 2, sf::Lines);
		ioWindow.draw(p.shape);

		//Pole debug
		/*
		float teta = 1.5;
		sf::Vector2f aRVect = unitVectorFromAngle(p.orientation+teta/2.0+M_PI/2.0);
		sf::Vector2f aRPole = p.position + 5.0f * aRVect;
		sf::Vector2f aLVect = unitVectorFromAngle(p.orientation-teta/2.0-M_PI/2.0);
		sf::Vector2f aLPole = p.position + 5.0f * aLVect;
		sf::CircleShape aRPoleDot(1);
		aRPoleDot.setOrigin(1, 1);
		aRPoleDot.setPosition(aRPole);
		sf::CircleShape aLPoleDot(1);
		aLPoleDot.setOrigin(1, 1);
		aLPoleDot.setPosition(aLPole);
		aRPoleDot.setFillColor(sf::Color::Red);
		aLPoleDot.setFillColor(sf::Color::Red);
		ioWindow.draw(aRPoleDot);
		ioWindow.draw(aLPoleDot);
		*/
	}
}

// Main function
int main()
{
	// Test
	//std::cout << angleFromVector(sf::Vector2f(-1.0,-2.0)) << std::endl;

	// Create the main window
	sf::RenderWindow window(sf::VideoMode(WORLD_WIDTH, WORLD_HEIGTH), "Particle system");
	ImGui::SFML::Init(window);

	// Model
	Model myModel;

	// Create a view with the same size as the window
	sf::View view(sf::FloatRect(-( ( float )WORLD_WIDTH ) * 0.5f, -( ( float )WORLD_HEIGTH ) * 0.5f, WORLD_WIDTH, WORLD_HEIGTH));

	// Variables to store the state of mouse dragging
	bool isDragging = false;
	sf::Vector2f lastMousePos;

	sf::Clock deltaClock;

	// main loop
	while ( window.isOpen() ) {
		// Handle events
		sf::Event event;
		while ( window.pollEvent(event) ) {
			ImGui::SFML::ProcessEvent(window, event);

			switch ( event.type )
			{
			case sf::Event::Closed:
				window.close();
				break;
			case sf::Event::MouseWheelScrolled:
				// If the wheel scrolled up, zoom in, else zoom out
				if ( event.mouseWheelScroll.delta > 0 )
					view.zoom(0.9f);
				else
					view.zoom(1.1f);
				break;
			case sf::Event::MouseButtonPressed:
				if ( event.mouseButton.button == sf::Mouse::Left ) {
					isDragging = true;
					lastMousePos = window.mapPixelToCoords(sf::Mouse::getPosition(window));
				}
				break;
			case sf::Event::MouseButtonReleased:
				if ( event.mouseButton.button == sf::Mouse::Left )
					isDragging = false;
				break;
			}
		}

		ImGui::SFML::Update(window, deltaClock.restart());

		ImGui::Begin("Hello, world!");
		if (ImGui::Button("Reset") )
		{
			myModel.init();
		}
		ImGui::DragFloat("Sigma", &sigmaAttraction, 0.1f, 0.001f, DOT_SIZEf, "%.5f", ImGuiSliderFlags_Logarithmic);
		ImGui::End();

		// Check for mouse dragging
		if ( isDragging && sf::Mouse::isButtonPressed(sf::Mouse::Left) ) {
			sf::Vector2f const mousePos = window.mapPixelToCoords(sf::Mouse::getPosition(window));
			sf::Vector2f const delta = lastMousePos - mousePos;

			view.move(delta);
		}

		// Refresh
		window.setView(view);

		// Update the last mouse position
		lastMousePos = window.mapPixelToCoords(sf::Mouse::getPosition(window));

		// Model update
		myModel.step();

		// Clear screen
		window.clear();

		// Draw particles
		drawModel(window, myModel);

		ImGui::SFML::Render(window);
		// Update the window
		window.display();
	}

	ImGui::SFML::Shutdown();

	return 0;
}
