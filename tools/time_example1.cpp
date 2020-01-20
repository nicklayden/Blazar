mp_type time_example1(mp_type z, mp_type eta)
{
	mp_type t1 , t2 , t3 , t8 , t9 , t11 , t12 , t13 , t16 , t23 ;


	t1 = pow(M, 0.1e1 / 0.3e1);
	t2 = t1 * t1;
	t3 = M * M;
	t8 = pow(0.1000000000e0 * t3 * M + 0.5000e4 * t3 + 0.1250000000e2, 0.1e1 / 0.3e1);
	t9 = t8 * t8;
	t11 = 0.1e1 / t9 * t2;
	t12 = sqrt(t11);
	t13 = t12 * t11;
	t16 = sin(eta);
	t23 = -0.1591549431e0 / t13 * (0.3141592654e5 * t13 * M + t16 - 0.1e1 * eta) * M;

	return t23 ;
}
