     for (unsigned int i=0; i<N-1; ++i) 
     {
         for (unsigned int j=i+2; j<N; ++j)               
         {
             t = 2*matr[i][j]/(matr[i][i] - matr[j][j]);
             phi = 0.5 * atan(t);
             c = cos(phi);
             s = sin(phi);

             bii = c*c*matr[i][i] + 2*c*s*matr[i][j] + s*s*matr[j][j];
             bij = s*c*(matr[j][j] - matr[i][i]) + matr[i][j] * (c*c - s*s);
             bjj = s*s*matr[i][i] + c*c*matr[j][j] - 2*c*s*matr[i][j];
             bji = bij;

             matr[i][i] = bii;
             matr[i][j] = bij;
             matr[j][i] = bji;
             matr[j][j] = bjj;
        }
    }
