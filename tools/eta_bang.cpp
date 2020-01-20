mp_type eta_bang(mp_type z, mp_type eta)
{
	mp_type t1 , t2 , t8 , t20 ;


	t1 = z * z;
	t2 = t1 * z;
	t8 = sin(eta);
	t20 = 0.1e1 / (0.1e0 * t2 + 0.5000e4 * t1 + 0.125e2) * (0.1e0 * t2 * eta + 0.5000e4 * t1 * eta + 0.125e2 * eta - 0.1e0 * t2 * t8 - 0.5000e4 * t1 * t8 - 0.125e2 * t8 - 0.3141592656e5 * t1);

	return t20 ;
}
