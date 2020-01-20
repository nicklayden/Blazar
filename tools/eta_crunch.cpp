mp_type eta_crunch(mp_type z, mp_type eta)
{
	mp_type t1 , t2 , t8 , t21 ;


	t1 = z * z;
	t2 = t1 * z;
	t8 = sin(eta);
	t21 = 0.1e1 / (0.1e0 * t2 + 0.5000e4 * t1 + 0.125e2) * (0.1e0 * t2 * eta + 0.5000e4 * t1 * eta + 0.125e2 * eta - 0.1e0 * t2 * t8 - 0.5000e4 * t1 * t8 - 0.125e2 * t8 - 0.6283185312e1 * t + 0.6283185312e0 * t2 + 0.785398164e2);

	return t21 ;
}
