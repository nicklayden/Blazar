mp_type app_horizon(mp_type z, mp_type eta)
{
	mp_type t4 , t5 , t6 , t10 , t11 , t14 , t19 ;


	t4 = z * z;
	t5 = t4 * t4;
	t6 = t5 * t5;
	t10 = (int) pow((double) (-2 * t6 * t4 + t6 + 1), (double) 2);
	t11 = t10 * t10;
	t14 = cos(eta);
	t19 = 0.3200000000e0 * (1 - t14) / t11 * z - t4 * z;

	return t19 ;
}
