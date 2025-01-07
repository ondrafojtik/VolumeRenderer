#include "stdafx.h"

#include "cell.h"
#include "matrix4x4.h"

Cell::Cell( const DATA_TYPE * rhos, const Vector3 & A, const Vector3 & B )
{
	memcpy( rhos_, rhos, sizeof( rhos_ ) );

	bounds_ = AABB( A, B );
}

Cell::~Cell()
{

}

float Cell::Gamma( const Vector3 & uvw ) const
{	
	//return ( rho_A() + rho_B() + rho_C() + rho_D() + rho_E() + rho_F() + rho_G() + rho_H() ) / 8.0f;
	//return tricubicInterpolate(rhos_, uvw.x, uvw.y, uvw.z);

	
	const float & u_t = uvw.x;
	const float & v_t = uvw.y;
	const float & w_t = uvw.z;
	
	// trilinear interpolation
	const COMP_TYPE alpha_AB = rho_A() * ( 1 - u_t ) + rho_B() * u_t;
	const COMP_TYPE alpha_DC = rho_D() * ( 1 - u_t ) + rho_C() * u_t;
	const COMP_TYPE alpha_EF = rho_E() * ( 1 - u_t ) + rho_F() * u_t;
	const COMP_TYPE alpha_HG = rho_H() * ( 1 - u_t ) + rho_G() * u_t;
	
	const COMP_TYPE beta_0 = alpha_AB * ( 1 - v_t ) + alpha_DC * v_t;
	const COMP_TYPE beta_1 = alpha_EF * ( 1 - v_t ) + alpha_HG * v_t;
	
	return static_cast< DATA_TYPE >( beta_0 * ( 1 - w_t ) + beta_1 * w_t );

	// TODO try to add tricubic interpolation
}



Matrix4x4 hx_{
	-1, 0, 1, 0,
	-2, 4, 2, 0,
	-1, 0, 1, 0,
	 0, 0, 0, 0
};

Matrix4x4 hy_{
	-1, -2, -1, 0,
	 0,  4,  0, 0,
	 1,  2,  1, 0,
	 0,  0,  0, 0
};

Matrix4x4 hz_{
	-1, -2, -1, 0,
	 0,  4,  0, 0,
	 1,  2,  1, 0,
	 0,  0,  0, 0
};


// https://dsp.stackexchange.com/questions/45433/how-sobel-edge-detector-as-a-first-order-derivative-is-turned-to-a-2-dim-filter
// https://en.wikipedia.org/wiki/Sobel_operator
Matrix4x4 kernel{
	1, 2, 1, 0,
	2, 4, 2, 0,
	1, 2, 1, 0,
	0, 0, 0, 0
};


Vector3 Cell::GradGamma(const Vector3& uvw) const
{
	// TODO compute the gradient of the scalar field here (use finite central differences or Sobel–Feldman operator)
	
	
	float Gx = 0;
	float Gy = 0;
	float Gz = 0;
	
	//float Gx_ = 0;
	//float Gy_ = 0;
	//float Gz_ = 0;


	for (int x = -1; x < 2; x++) 
		for (int y = -1; y < 2; y++) 
			for (int z = -1; z < 2; z++)
			{
				Vector3 p = uvw + Vector3(x, y, z);
				float val = Gamma(p);

				Gx += x * kernel.get(y, z) * val;
				Gy += y * kernel.get(x, z) * val;
				Gz += z * kernel.get(x, y) * val;

				//Gx_ += x * hx_.get(y, z) * val;
				//Gy_ += y * hy_.get(x, z) * val;
				//Gz_ += z * hz_.get(x, y) * val;

			}
			
	Vector3 result_(Gx, Gy, Gz);
	return result_;
	
}

float Cell::Integrate( Ray & ray, const float t0, const float t1 ) const
{	
	// TODO approximate computation of an integral of scalar values along the given segment of the ray
	
	float values = 0.0f;

	float step = 0.1f;
	for (float i = t0; i <= t1; i += step)
	{
		float value_ = Gamma(u(ray.eval(i)));
		values += value_ * step; // value * dt ?
	}
	
	return values;

}

float Cell::FindIsoSurface( Ray & ray, const float t0, const float t1, const float iso_value ) const
{	
	// TODO find the parametric distance of the iso surface of the certain iso value along the given segment of the ray
	float distance_ = -1.0f;
	float step = 0.1f;
	
	for (float i = t0; i <= t1; i += step)
	{
		float val = Gamma(u(ray.eval(i)));

		if (val >= iso_value) 
		{
			//distance_ = val;
			distance_ = i;
			break;
		}
	}
	return distance_;
	
}

float Cell::FrontToBack(Ray& ray, float t0, float t1, Vector3* current_color, float* sample_alpha)
{

	float step = 0.01f;
	float is_1 = 0;

	float c_out = current_color->x;
	float a_out = *sample_alpha;

	float skin_th = 0.3f;
	float bone_th = 0.50f;

	for (float i = t0; i <= t1; i += step)
	{
		float C = Gamma(u(ray.eval(i)));
		if (C < skin_th) continue;

		//skin
		float a = 0.00055f;
		if (C > bone_th) a = 0.01f;
		
		c_out += ((1 - a_out) * a * C);
		a_out += ((1 - a_out) * a);

		if (a_out >= 1.0f)
		{
			is_1 = 1;
			break;
		}
	}

	*current_color = Vector3(c_out, c_out, c_out);
	*sample_alpha = a_out;

	return is_1;
	
}

float Cell::FrontToBack2(Ray& ray, float t0, float t1, Vector3* current_color, float* sample_alpha)
{
	float step = 0.1f;
	float is_1 = 0;

	float c_out = current_color->x;
	float a_out = *sample_alpha;

	float skin_th = 0.3f;
	float bone_th = 0.50f;

	Vector3 ccin = *current_color;
	float ccinx = current_color->x;

	for (float i = t0; i <= t1; i += step)
	{
		float C = Gamma(u(ray.eval(i)));
		if (C < skin_th) continue;

		//skin
		float a = 0.0015f;
		if (C > bone_th) a = 0.3f;

		c_out += ((1 - a_out) * a * C);
		a_out += ((1 - a_out) * a);

		auto cc = GradGamma(u(ray.eval(i)));
		ccin = ccin + ((1 - a_out) * a * cc);
		ccinx = ccinx + ((1 - a_out) * a * ccin.x);

		if (a_out >= 1.0f)
		{
			is_1 = 1;
			break;
		}
	}
	//float final__ = (0.2989 * ggin.x) + (0.587 * ggin.y) + (0.114 * ggin.z);
	float final__ = (0.33 * ccin.x) + (0.33 * ccin.y) + (0.33 * ccin.z);


	//*current_color = Vector3(gginx, gginx, gginx);
	//*current_color = ggin;
	*current_color = Vector3(final__, final__, final__);
	*sample_alpha = a_out;

	return is_1;

}

float Cell::BackToFront(Ray& ray, float t0, float t1, Vector3* current_color, float* sample_alpha)
{

	float step = 0.01f;
	float is_1 = 0;

	float c_out = current_color->x;
	float a_out = *sample_alpha;

	float skin_th = 0.3f;
	float bone_th = 0.50f;

	for (float i = t0; i <= t1; i += step)
	{
		float C = Gamma(u(ray.eval(i)));
		if (C < skin_th) continue;

		//skin
		float a = 0.00055f;
		if (C > bone_th) a = 0.01f;


		c_out = ((1 - a) * c_out) + (a * C);
		a_out = ((1 - a) * a_out) + a;

		if (a_out >= 1.0f)
		{
			is_1 = 1;
			break;
		}
	}

	*current_color = Vector3(c_out, c_out, c_out);
	*sample_alpha = a_out;

	return is_1;

}

Vector3 Cell::u( const Vector3 & p ) const
{
	Vector3 uvw = ( p - A() ) / ( G() - A() ); // gives the reference coordinates of the world space point p inside this cell (uvw is in the range <0,1>^3)

	return uvw;
}

Vector3 Cell::A() const
{
	return bounds_.lower_bound();
}

Vector3 Cell::G() const
{
	return bounds_.upper_bound();
}

float Cell::rho_A() const
{
	return rhos_[0];
}

float Cell::rho_B() const
{
	return rhos_[1];
}

float Cell::rho_C() const
{
	return rhos_[2];
}

float Cell::rho_D() const
{
	return rhos_[3];
}

float Cell::rho_E() const
{
	return rhos_[4];
}

float Cell::rho_F() const
{
	return rhos_[5];
}

float Cell::rho_G() const
{
	return rhos_[6];
}

float Cell::rho_H() const
{
	return rhos_[7];
}
