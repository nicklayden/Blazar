mp_type eta_eq(mp_type z, mp_type eta)
{
	mp_type t1 , t6 , t7 , t9 , t11 , t12 , t14 , t16 , t20 ;


	t1 = z * z;
	t6 = pow(0.1e0 * t1 * z + 0.5000e4 * t1 + 0.125e2, 0.1e1 / 0.3e1);
	t7 = t6 * t6;
	t9 = sin(eta);
	t11 = pow(z, 0.1e1 / 0.3e1);
	t12 = t11 * t11;
	t14 = 0.1e1 / t7;
	t16 = sqrt(t14 * t12);
	t20 = t14 * (t7 * eta - t7 * t9 - 0.3141592656e5 * t16 * t12 * z);

	return t20 ;
}
