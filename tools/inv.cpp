mp_50 inv(mp_50 z, mp_50 eta)
{
	mp_50 t4 , t5 , t6 , t9 , t10 , t11 , t13 , t14 , t21 , t23 , t26 , t27 , t30 , t35 , t36 , t37 , t38 , t40 , t41 , t43 ;
	mp_50 t46 , t52 , t69 ;


	t4 = (int) ((double) z * (double) z);
	t5 = (int) ((double) t4 * (double) t4);
	t6 = (int) ((double) t5 * (double) t5);
	t9 = -2 * t6 * t4 + t6 + 1;
	t10 = t9 * t9;
	t11 = t10 * t10;
	t13 = cos(eta);
	t14 = 0.1e1 - t13;
	t21 = sqrt(t13 * t14);
	t23 = sqrt((double) t11);
	t26 = 0.1000000000e1 * t13;
	t27 = acos(t26);
	t30 = sqrt((0.1000000000e1 + t26) * t14);
	t35 = abs(-t9);
	t36 = t35 * t35;
	t37 = t36 * t36;
	t38 = sqrt((double) t37);
	t40 = 0.1e1 / t38 / (double) t37;
	t41 = 0.6250000000e1 * (double) t37;
	t43 = acos(-0.1e1 + t41);
	t46 = sqrt(0.1250000000e2 - 0.3906250000e2 * (double) t37);
	t52 = sqrt(0.2e1 - t41);
	t69 = 0.6400000000e0 * t14 / (double) t11 + 0.1060660172e2 * (double) (8 * t5 * t4 * z - 20 * t6 * z) * (0.2560000000e0 * (t27 - 0.9999999999e0 * t30) / t23 / (double) t11 + 0.9050966799e-1 * (-t46 * (double) t36 + 0.3141592654e1 - t43) * t40 - 0.9050966799e-1 * (0.3141592654e1 + t43 - 0.2500000000e1 * t52 * (double) t36) * t40) * t23 * t21 / t14 / (double) t9;
	return t69 ;
}
