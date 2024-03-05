Omega_AI +=  +0.50000000  * extract_mat(h_p, "AI", o, v);
Omega_AI = Omega_AI .+  -1.00000000  * np.einsum("ii->", extract_mat(g_p, "AIoo", o, v), optimize="optimal");
Omega_AI = Omega_AI .+  +1.00000000  * np.einsum("ia,ai->", extract_mat(F, "ov", o, v), extract_mat(s, "AIvo", o, v), optimize="optimal");
Omega_AI = Omega_AI .+  -1.00000000  * np.einsum("Bia,Bai->", extract_mat(g_p, "AVov", o, v), extract_mat(s, "VIvo", o, v), optimize="optimal");
Omega_AI = Omega_AI .+  +0.50000000  * np.einsum("iajb,aibj->", extract_mat(L, "ovov", o, v), extract_mat(s2, "AIvovo", o, v), optimize="optimal");
