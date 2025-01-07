#include "stdafx.h"

#include "camera.h"

Ray Camera::GenerateRay( const float sx, const float sy )
{		
	const float x = pixel_size_ * ( sx - static_cast<float>( 0.5 ) * ( width_ - 1 ) );
	const float y = pixel_size_ * ( -sy + static_cast<float>( 0.5 ) * ( height_ - 1 ) );

	Vector3 direction = Vector3( x, y, -d_ ); // sm�r nov�ho paprsku v kamerov�m prostoru	
	direction.Normalize();

	Vector3 direction_transformed = Vector3(
		direction.x * axis_x_.x + direction.y * axis_y_.x + direction.z * axis_z_.x,
		direction.x * axis_x_.y + direction.y * axis_y_.y + direction.z * axis_z_.y,
		direction.x * axis_x_.z + direction.y * axis_y_.z + direction.z * axis_z_.z
		); // sm�r nov�ho paprsku ve sv�tov�m prostoru
	direction_transformed.Normalize();

	return Ray( view_from_, direction_transformed, 0, 0 );
}
