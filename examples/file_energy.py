E +=  +1.00000000  * extract_mat(h_p, "II", o, v);
E = E .+  +1.00000000  * np.einsum("ii->", extract_mat(F, "oo", o, v), optimize="optimal");
E = E .+  -2.00000000  * np.einsum("ii->", extract_mat(g_p, "IIoo", o, v), optimize="optimal");
E = E .+  +1.00000000  * np.einsum("ii->", extract_mat(h, "oo", o, v), optimize="optimal");
