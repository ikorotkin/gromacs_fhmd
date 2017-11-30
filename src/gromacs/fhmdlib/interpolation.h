#ifndef FHMD_INTERPOLATION_H_
#define FHMD_INTERPOLATION_H_

void trilinear_find_neighbours(const rvec x, const int n, dvec xi, int *nbr, FHMD *fh);
void trilinear_interpolation(dvec result, dvec xi, dvec f0, dvec f1, dvec f2, dvec f3, dvec f4, dvec f5, dvec f6, dvec f7);

#define INTERPOLATE(f) arr[nbr[0]].f, arr[nbr[1]].f, arr[nbr[2]].f, arr[nbr[3]].f, arr[nbr[4]].f, arr[nbr[5]].f, arr[nbr[6]].f, arr[nbr[7]].f

#endif /* FHMD_INTERPOLATION_H_ */
