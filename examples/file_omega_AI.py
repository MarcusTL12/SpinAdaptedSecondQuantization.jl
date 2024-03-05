Omega_AI_ +=  +0.50000000  * extract_mat(h_p, "ii", o, v);
Omega_AI_ = Omega_AI_ .+  -1.00000000  * np.einsum("ii->", extract_mat(g_p, "iioo", o, v), optimize="optimal");
Omega_AI_ = Omega_AI_ .+  +1.00000000  * np.einsum("ia,ai->", extract_mat(F, "oo", o, v), extract_mat(s, "iioo", o, v), optimize="optimal");
Omega_AI_ = Omega_AI_ .+  -1.00000000  * np.einsum("Bia,Bai->", extract_mat(g_p, "iooo", o, v), extract_mat(s, "oioo", o, v), optimize="optimal");
Omega_AI_ = Omega_AI_ .+  +0.50000000  * np.einsum("iajb,aibj->", extract_mat(L, "oooo", o, v), extract_mat(s, "iioooo", o, v), optimize="optimal");
