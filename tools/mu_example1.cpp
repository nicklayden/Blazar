mp_type mu_example1(mp_type M, mp_type eta)
{
	mp_type t1 , t2 , t3 , t5 , t6 , t8 , t13 , t14 , t15 , t17 , t18 , t20 , t26 , t28 , t37 ;


	t1 = pow(2, 0.1e1 / 0.6e1);
	t2 = pow(2, 0.1e1 / 0.3e1);
	t3 = t2 * t2;
	t5 = pow(0.3141592654e1 * M, 0.1e1 / 0.3e1);
	t6 = t5 * t5;
	t8 = M * M;
	t13 = pow(0.1e0 * t8 * M + 5000 * t8 + 0.125e2, 0.1e1 / 0.3e1);
	t14 = t13 * t13;
	t15 = 0.1e1 / t14;
	t17 = 0.10e1 * t15 * t6 * t3;
	t18 = cos(eta);
	t20 = 0.1e1 / (1 - t18);
	t26 = sqrt(-t17 + 0.2000000000e1 * t15 * t6 * t3 * t20);
	t28 = sqrt(1 - t17);
	t37 = -0.1000000000e1 * t15 * t6 * t20 / M * (t26 + t28) * t1;

	return t37 ;
}
