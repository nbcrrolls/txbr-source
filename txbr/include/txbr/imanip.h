#ifndef IMANIP_H
#define IMANIP_H 

void translate_image(int nx, int ny, float* data_2D, int t_x, int t_y);
void rotate_image(int nx, int ny, float* data_2D, double center_x, double center_y, double theta);

#endif
