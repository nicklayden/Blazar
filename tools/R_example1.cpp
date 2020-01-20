mp_type R_example1(mp_type z, mp_type eta)
{
	mp_type t1 , t4 , t6 , t7 , t10 , t15 , t16 , t21 ;


	t1 = cos(eta);
	t4 = pow(2, 0.1e1 / 0.3e1);
	t6 = pow(0.3141592654e1 * M, 0.1e1 / 0.3e1);
	t7 = t6 * t6;
	t10 = M * M;
	t15 = pow(0.1e0 * t10 * M + 5000 * t10 + 0.125e2, 0.1e1 / 0.3e1);
	t16 = t15 * t15;
	t21 = 0.5000000000e0 * t16 / t7 * t4 * (1 - t1) * M - 2 * M;

	return t21 ;
}
