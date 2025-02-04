#include "point.hpp"

xy_tuple getxy(CodiPoint self)
{
	xy_tuple xy = std::make_tuple(self.x, self.y);
	return xy;
}

void setNormals(CodiPoint* globaldata, int idx, xy_tuple n)
{
	globaldata[idx].nx = std::get<0>(n);
	globaldata[idx].ny = std::get<1>(n);
}