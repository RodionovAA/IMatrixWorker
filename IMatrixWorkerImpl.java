import java.util.Arrays;

interface IMatrixWorker {
   public void print(double[][] m);
   public boolean haveSameDimension(double[][] m1, double[][] m2);
   public  double[][] add(double[][] m1, double[][] m2);
   public double[][] subtract(double[][] m1, double[][] m2);
   public  double[][] multiply(double[][] m1, double[][] m2);
}

public class IMatrixWorkerImpl implements IMatrixWorker  {
    
    private double[][] matrix;
    
    public void setMatrix(double[][] m){
        matrix = new double[m.length][m[0].length];
        for(int i=0;i<m.length;i++){
           for(int j=0;j<m[0].length;j++){
              matrix[i][j] = m[i][j]; 
           }     
        }
    }
    
    public double[][] getMatrix(){
        return matrix;
    }
    
    public IMatrixWorkerImpl(int N, int M){
        matrix = new double[N][M];
        for(int i=0;i<N;i++){
           for(int j=0;j<M;j++){
              matrix[i][j] = 0.0; 
           }     
        }
    }
    
    @Override
    public void print(double[][] m) {
        for (int i = 0;i< m.length;i++){
            for (int j = 0; j<m[i].length;j++){
                 if (j == m[i].length - 1){
                  System.out.println(Double.toString(m[i][0])+ " ");   
                } else {
                System.out.print(Double.toString(m[i][0])+ " ");
                }
            }
        }
    }

    @Override
    public boolean haveSameDimension(double[][] m1, double[][] m2) {
        boolean res = false;
        //исхожу из того, что матриц с разной размерностью столбцов быть не может
        if ((m1.length == m2.length) && (m1[0].length == m2[0].length)){
            res = true;
        }
        return res;
    }

    @Override
    public  double[][] add(double[][] m1, double[][] m2) {
        double[][] res = new double[m1.length][];
        for(int i = 0; i < res.length; i++) {
            res[i] = Arrays.copyOf(m1[i],m1.length);
            for (int j = 0; j < res[i].length; j++) {
                res[i][j] += m2[i][j];
            }
        }
        return res;
    }

    @Override
    public double[][] subtract(double[][] m1, double[][] m2) {
        double[][] res = new double[m1.length][];
        for(int i = 0; i < res.length; i++) {
            res[i] = Arrays.copyOf(m1[i],m1.length);
            for (int j = 0; j < res[i].length; j++) {
                res[i][j] -= m2[i][j];
            }
        }
        return res;
    }

    @Override
    public  double[][] multiply(double[][] m1, double[][] m2) {
        if( m1[0].length != m2.length){
            throw new IllegalArgumentException("Sorry, input data is incorrect");
        }
        double[][] res = new double[m1.length][m2[0].length];
        
        
        for(int i = 0; i< res.length; i++){
            for(int j = 0;j<res[0].length;j++){
                res[i][j] = 0;
                for(int k = 0; k < m1.length;k++){
                    res[i][j] += m1[i][k] * m2[k][j];
                }
            }
        }
        return res;
    }
    
    private double[][] generateSubArray (double A[][], int N, int j1){
        double m[][]; 
        m = new double[N-1][];
            for (int k=0; k<(N-1); k++)
                    m[k] = new double[N-1];

            for (int i=1; i<N; i++){
                  int j2=0;
                  for (int j=0; j<N; j++){
                      if(j == j1)
                            continue;
                      m[i-1][j2] = A[i][j];
                      j2++;
                  }
              }
            return m;
    }

    public double determinant(double A[][], int N){
        if(A.length != A[0].length){
           //определитель вычисляем только для квадратных матриц
           // опять же исхожу из того, что столбцов с разной размерностью быть не может
            throw new IllegalArgumentException("Sorry, input array is incorrect");
        }
        double res;
        double m[][];

        if (N == 1){
            res = A[0][0];
        } else if (N == 2) {
            res = A[0][0]*A[1][1] - A[1][0]*A[0][1];
        } else{
            res=0;
            for (int j1=0; j1<N; j1++){
                 m = generateSubArray (A, N, j1);
                 res += Math.pow(-1.0, 1.0+j1+1.0) * A[0][j1] * determinant(m, N-1);
            }
        }
        return res;
    } 
}
