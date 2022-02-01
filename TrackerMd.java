
import java.io.File;
import java.io.FileNotFoundException;
import java.util.Locale;
import java.util.Scanner;

public class TrackerMd {
    
    // initial values
    private Matrix A;   // velocity transformation matrix 
    private Matrix B;   // acceleration transformation matrix
    private Matrix stateMatrix;   // X 
    private Matrix controlVariableMatrix;    // u
    private Matrix noiseInProcessMatrix;     // w
    private Matrix errorInCalculatingTheProcessCovarianceMatrix;  // Q
    private Matrix H;   // for matrix shape compatibility
    private Matrix processCovarianceMatrix;  // P
    private Matrix observationCovarianceMatrix; // R
    private Matrix C;   // for observation transformation
    private Matrix Z;   // errors due to the sensor manifacturing electronics goes here
    private Matrix identityMatrix;  // I
    
    // dynamic values
    private Matrix predictedStateMatrix;
    private Matrix predictedProcessCovarianceMatrix; // P
    private Matrix kalmanGain;         //K
    private Matrix observationMatrix;  // Y
    private Matrix currentState;    // X
    private Matrix updatedProcessCovarianceMatrix; // P
    private Matrix observationStateMatrix;  // Ym
    
    TrackerMd(double initialPositionX,     double initialPositionY,    double initialPositionZ,
              double intitialVelocityX,    double intitialVelocityY,    double intitialVelocityZ,
              double initialAccelerationX, double initialAccelerationY, double initialAccelerationZ, 
              double deltaT, 
              double processErrorInX,      double processErrorInY,      double processErrorInZ,
              double processErrorInVx,     double processErrorInVy,     double processErrorInVz,
              double observationErrorInX,  double observationErrorInY,  double observationErrorInZ, 
              double observationErrorInVx, double observationErrorInVy, double observationErrorInVz){
        
        A = new Matrix (new double[][]{
            {1.0, 0.0, 0.0, deltaT, 0.0   , 0.0   },
            {0.0, 1.0, 0.0, 0.0   , deltaT, 0.0   },
            {0.0, 0.0, 1.0, 0.0   , 0.0   , deltaT},
            {0.0, 0.0, 0.0, 1.0   , 0.0   , 0.0   },
            {0.0, 0.0, 0.0, 0.0   , 1.0   , 0.0   },
            {0.0, 0.0, 0.0, 0.0   , 0.0   , 1.0   }
        });
        B = new Matrix (new double[][]{
            {deltaT * deltaT /2.0, 0.0                 , 0.0                 },
            {0.0                 , deltaT * deltaT /2.0, 0.0                 },
            {0.0                 , 0.0                 , deltaT * deltaT /2.0},
            {deltaT              , 0.0                 , 0.0                 },
            {0.0                 , deltaT              , 0.0                 },
            {0.0                 , 0.0                 , deltaT              }
        });
        stateMatrix=new Matrix (new double[][]{
            {initialPositionX},
            {initialPositionY},
            {initialPositionZ},
            {intitialVelocityX},
            {intitialVelocityY},
            {intitialVelocityZ}
        });
        controlVariableMatrix = new Matrix (new double[][]{
            {initialAccelerationX},
            {initialAccelerationY},
            {initialAccelerationZ}
        });
        noiseInProcessMatrix = new Matrix (new double[][]{
            {0.0},
            {0.0},
            {0.0},
            {0.0},
            {0.0},
            {0.0}
        });

        errorInCalculatingTheProcessCovarianceMatrix = new Matrix(new double[][]{
            {0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
        });     // Q matrix   zero for simplicity

        H= new Matrix(new double[][]{
            {1.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {0.0, 1.0, 0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 1.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 1.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0, 1.0, 0.0},
            {0.0, 0.0, 0.0, 0.0, 0.0, 1.0}
        }); // identity matrix 

        C = new Matrix (new double[][]{
            {1.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {0.0, 1.0, 0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 1.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 1.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0, 1.0, 0.0},
            {0.0, 0.0, 0.0, 0.0, 0.0, 1.0}
        }); // identity matrix

        Z = new Matrix(new double[][]{  // errors due to the sensor manifacturing electronics goes here
            {0.0},
            {0.0},
            {0.0},
            {0.0},
            {0.0},
            {0.0}
        });

        identityMatrix = new Matrix(new double[][]{
            {1.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {0.0, 1.0, 0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 1.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 1.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0, 1.0, 0.0},
            {0.0, 0.0, 0.0, 0.0, 0.0, 1.0}
        }); 
        
        processCovarianceMatrix = new Matrix(new double[][]{{processErrorInX*processErrorInX ,processErrorInY*processErrorInX ,processErrorInZ*processErrorInX  ,processErrorInVx*processErrorInX  ,processErrorInVy*processErrorInX  ,processErrorInVz*processErrorInX },
                                                            {processErrorInX*processErrorInY ,processErrorInY*processErrorInY ,processErrorInZ*processErrorInY  ,processErrorInVx*processErrorInY  ,processErrorInVy*processErrorInY  ,processErrorInVz*processErrorInY },
                                                            {processErrorInX*processErrorInZ ,processErrorInY*processErrorInZ ,processErrorInZ*processErrorInZ  ,processErrorInVx*processErrorInZ  ,processErrorInVy*processErrorInZ  ,processErrorInVz*processErrorInZ },
                                                            {processErrorInX*processErrorInVx,processErrorInY*processErrorInVx,processErrorInZ*processErrorInVx ,processErrorInVx*processErrorInVx ,processErrorInVy*processErrorInVx ,processErrorInVz*processErrorInVx},
                                                            {processErrorInX*processErrorInVy,processErrorInY*processErrorInVy,processErrorInZ*processErrorInVy ,processErrorInVx*processErrorInVy ,processErrorInVy*processErrorInVy ,processErrorInVz*processErrorInVy},
                                                            {processErrorInX*processErrorInVz,processErrorInY*processErrorInVz,processErrorInZ*processErrorInVz ,processErrorInVx*processErrorInVz ,processErrorInVy*processErrorInVz ,processErrorInVz*processErrorInVz}});      
        processCovarianceMatrix.setValueAt(0,1,0.0);    // for simplicity
        processCovarianceMatrix.setValueAt(0,2,0.0);    // for simplicity
        processCovarianceMatrix.setValueAt(0,3,0.0);    // for simplicity
        processCovarianceMatrix.setValueAt(0,4,0.0);    // for simplicity
        processCovarianceMatrix.setValueAt(0,5,0.0);    // for simplicity
        
        processCovarianceMatrix.setValueAt(1,0,0.0);    // for simplicity
        processCovarianceMatrix.setValueAt(1,2,0.0);    // for simplicity
        processCovarianceMatrix.setValueAt(1,3,0.0);    // for simplicity
        processCovarianceMatrix.setValueAt(1,4,0.0);    // for simplicity
        processCovarianceMatrix.setValueAt(1,5,0.0);    // for simplicity
        
        processCovarianceMatrix.setValueAt(2,0,0.0);    // for simplicity
        processCovarianceMatrix.setValueAt(2,1,0.0);    // for simplicity
        processCovarianceMatrix.setValueAt(2,3,0.0);    // for simplicity
        processCovarianceMatrix.setValueAt(2,4,0.0);    // for simplicity
        processCovarianceMatrix.setValueAt(2,5,0.0);    // for simplicity
        
        processCovarianceMatrix.setValueAt(3,0,0.0);    // for simplicity
        processCovarianceMatrix.setValueAt(3,1,0.0);    // for simplicity
        processCovarianceMatrix.setValueAt(3,2,0.0);    // for simplicity
        processCovarianceMatrix.setValueAt(3,4,0.0);    // for simplicity
        processCovarianceMatrix.setValueAt(3,5,0.0);    // for simplicity
        
        processCovarianceMatrix.setValueAt(4,0,0.0);    // for simplicity
        processCovarianceMatrix.setValueAt(4,1,0.0);    // for simplicity
        processCovarianceMatrix.setValueAt(4,2,0.0);    // for simplicity
        processCovarianceMatrix.setValueAt(4,3,0.0);    // for simplicity
        processCovarianceMatrix.setValueAt(4,5,0.0);    // for simplicity
        
        processCovarianceMatrix.setValueAt(5,0,0.0);    // for simplicity
        processCovarianceMatrix.setValueAt(5,1,0.0);    // for simplicity
        processCovarianceMatrix.setValueAt(5,2,0.0);    // for simplicity
        processCovarianceMatrix.setValueAt(5,3,0.0);    // for simplicity
        processCovarianceMatrix.setValueAt(5,4,0.0);    // for simplicity
        
        observationCovarianceMatrix = new Matrix(new double[][]{{observationErrorInX*observationErrorInX ,observationErrorInY*observationErrorInX ,observationErrorInZ*observationErrorInX  ,observationErrorInVx*observationErrorInX  ,observationErrorInVy*observationErrorInX  ,observationErrorInVz*observationErrorInX },
                                                                {observationErrorInX*observationErrorInY ,observationErrorInY*observationErrorInY ,observationErrorInZ*observationErrorInY  ,observationErrorInVx*observationErrorInY  ,observationErrorInVy*observationErrorInY  ,observationErrorInVz*observationErrorInY },
                                                                {observationErrorInX*observationErrorInZ ,observationErrorInY*observationErrorInZ ,observationErrorInZ*observationErrorInZ  ,observationErrorInVx*observationErrorInZ  ,observationErrorInVy*observationErrorInZ  ,observationErrorInVz*observationErrorInZ },
                                                                {observationErrorInX*observationErrorInVx,observationErrorInY*observationErrorInVx,observationErrorInZ*observationErrorInVx ,observationErrorInVx*observationErrorInVx ,observationErrorInVy*observationErrorInVx ,observationErrorInVz*observationErrorInVx},
                                                                {observationErrorInX*observationErrorInVy,observationErrorInY*observationErrorInVy,observationErrorInZ*observationErrorInVy ,observationErrorInVx*observationErrorInVy ,observationErrorInVy*observationErrorInVy ,observationErrorInVz*observationErrorInVy},
                                                                {observationErrorInX*observationErrorInVz,observationErrorInY*observationErrorInVz,observationErrorInZ*observationErrorInVz ,observationErrorInVx*observationErrorInVz ,observationErrorInVy*observationErrorInVz ,observationErrorInVz*observationErrorInVz}});
        observationCovarianceMatrix.setValueAt(0,1,0.0);    // for simplicity
        observationCovarianceMatrix.setValueAt(0,2,0.0);    // for simplicity
        observationCovarianceMatrix.setValueAt(0,3,0.0);    // for simplicity
        observationCovarianceMatrix.setValueAt(0,4,0.0);    // for simplicity
        observationCovarianceMatrix.setValueAt(0,5,0.0);    // for simplicity
        
        observationCovarianceMatrix.setValueAt(1,0,0.0);    // for simplicity
        observationCovarianceMatrix.setValueAt(1,2,0.0);    // for simplicity
        observationCovarianceMatrix.setValueAt(1,3,0.0);    // for simplicity
        observationCovarianceMatrix.setValueAt(1,4,0.0);    // for simplicity
        observationCovarianceMatrix.setValueAt(1,5,0.0);    // for simplicity
        
        observationCovarianceMatrix.setValueAt(2,0,0.0);    // for simplicity
        observationCovarianceMatrix.setValueAt(2,1,0.0);    // for simplicity
        observationCovarianceMatrix.setValueAt(2,3,0.0);    // for simplicity
        observationCovarianceMatrix.setValueAt(2,4,0.0);    // for simplicity
        observationCovarianceMatrix.setValueAt(2,5,0.0);    // for simplicity
        
        observationCovarianceMatrix.setValueAt(3,0,0.0);    // for simplicity
        observationCovarianceMatrix.setValueAt(3,1,0.0);    // for simplicity
        observationCovarianceMatrix.setValueAt(3,2,0.0);    // for simplicity
        observationCovarianceMatrix.setValueAt(3,4,0.0);    // for simplicity
        observationCovarianceMatrix.setValueAt(3,5,0.0);    // for simplicity
        
        observationCovarianceMatrix.setValueAt(4,0,0.0);    // for simplicity
        observationCovarianceMatrix.setValueAt(4,1,0.0);    // for simplicity
        observationCovarianceMatrix.setValueAt(4,2,0.0);    // for simplicity
        observationCovarianceMatrix.setValueAt(4,3,0.0);    // for simplicity
        observationCovarianceMatrix.setValueAt(4,5,0.0);    // for simplicity
        
        observationCovarianceMatrix.setValueAt(5,0,0.0);    // for simplicity
        observationCovarianceMatrix.setValueAt(5,1,0.0);    // for simplicity
        observationCovarianceMatrix.setValueAt(5,2,0.0);    // for simplicity
        observationCovarianceMatrix.setValueAt(5,3,0.0);    // for simplicity
        observationCovarianceMatrix.setValueAt(5,4,0.0);    // for simplicity

    }

    private void findPredictedStateMatrix(){ // (step 1)     Xkp = A.Xk-1 + B.Uk + wk
        try{
            predictedStateMatrix = MatrixMathematics.add(
                                          MatrixMathematics.add(
                                                  MatrixMathematics.multiply(A, stateMatrix), 
                                                  MatrixMathematics.multiply(B, controlVariableMatrix)),
                                    noiseInProcessMatrix);
        }catch(IllegalDimensionException e){
             e.printStackTrace();
        }
    }
     
    
     private void findPredictedProcessCovarianceMatrix(){  // (step 3)   Pkp = A.Pk-1.AT + Qk
      try{
        predictedProcessCovarianceMatrix =  MatrixMathematics.add(
                                                 MatrixMathematics.multiply(
                                                         MatrixMathematics.multiply(A,processCovarianceMatrix),
                                                         MatrixMathematics.transpose(A)), 
                                              errorInCalculatingTheProcessCovarianceMatrix);
        
        predictedProcessCovarianceMatrix.setValueAt(0,1,0.0);    // for simplicity
        predictedProcessCovarianceMatrix.setValueAt(1,0,0.0);    // for simplicity
      }
      catch (IllegalDimensionException e){
        e.printStackTrace();
      }
    }
    
     
    private void calculateKalmanGain(){ // (step 4)       K = Pkp.HT / ( H.Pkp.HT + R )
       
       try{       
        kalmanGain = MatrixMathematics.multiply(MatrixMathematics.multiply(predictedProcessCovarianceMatrix, MatrixMathematics.transpose(H)),
                            MatrixMathematics.inverse(
                                    MatrixMathematics.add(
                                            MatrixMathematics.multiply(
                                                MatrixMathematics.multiply(H, predictedProcessCovarianceMatrix),
                                                MatrixMathematics.transpose(H)
                                            ),
                                             observationCovarianceMatrix
                                     )
                            ));
       }
        catch(IllegalDimensionException e){
            e.printStackTrace();
        }
        catch(NoSquareException e){
            e.printStackTrace();
        }
    }
    
    private void findObservationMatrix(){ // (step 5)     Yk = C.Ykm + Zm
            try{
            observationMatrix = MatrixMathematics.add(MatrixMathematics.multiply(C, observationStateMatrix), Z);
            }
            catch (IllegalDimensionException e){
             e.printStackTrace();
            }
    }
    
    
    private void calculateTheCurrentState(){ // (step 6)   Xk = Xkp + K[Yk-H.Xpk]    filtered state
        try{          
        currentState = MatrixMathematics.add(
                            predictedStateMatrix,
                            MatrixMathematics.multiply(
                                    kalmanGain,
                                    MatrixMathematics.subtract(
                                       observationMatrix,
                                       MatrixMathematics.multiply(H, predictedStateMatrix)
                                    )
                            )
                    );
        }
        catch(IllegalDimensionException e){
            e.printStackTrace();
        }
    }
    
    
    private void updateTheProcessCovarianceMatrix(){ // (step 7)  Pk = (I - K.H)Pkp
        try{
            updatedProcessCovarianceMatrix = MatrixMathematics.multiply(
                    MatrixMathematics.subtract(identityMatrix,
                            MatrixMathematics.multiply(kalmanGain, H)),
                            predictedProcessCovarianceMatrix
            );
        }catch(IllegalDimensionException e){
            e.printStackTrace();
        }
    }
    
    
    public void trackingStep(double ObservedX, double ObservedY, double ObservedZ, 
                             double ObservedVx, double ObservedVy, double ObservedVz){
        // step 1
        findPredictedStateMatrix();
        // step 3
        findPredictedProcessCovarianceMatrix();
        // step 4
        calculateKalmanGain();
        // step 5                  
        observationStateMatrix = new Matrix(new double[][]{
              {ObservedX},
              {ObservedY},
              {ObservedZ},
              {ObservedVx},
              {ObservedVy},
              {ObservedVz}
           });
        findObservationMatrix();        
        // step 6         
        calculateTheCurrentState();
        // step 7
        updateTheProcessCovarianceMatrix();
        // step 8  (current becomes previous)
        stateMatrix = currentState;
        processCovarianceMatrix = updatedProcessCovarianceMatrix;
    }
    
    
    public Matrix getPredictedStateMatrix(){
        return predictedStateMatrix;
    }
    public Matrix getPredictedProcessCovarianceMatrix(){
        return predictedProcessCovarianceMatrix;
    }
    public Matrix getKalmanGain(){
        return kalmanGain;
    }
    public Matrix getObservationMatrix(){
        return observationMatrix;
    }
    public Matrix getCurrentState(){
        return currentState;
    }
    public Matrix getUpdatedProcessCovarianceMatrix(){
        return updatedProcessCovarianceMatrix;
    }
    
    public static void main(String args[]) 
    {          
        double initialPositionX = 4000.0;     // X0   m
        double initialPositionY = 3000.0;     // Y0   m
        double initialPositionZ = 5000.0;     // Z0   m

        double intitialVelocityX = 280.0;     // V0x  m/s
        double intitialVelocityY = 250.0;     // V0y  m/s
        double intitialVelocityZ = 260.0;     // V0z  m/s

        double initialAccelerationX = 2.0;    // m/s2
        double initialAccelerationY = 3.0;    // m/s2
        double initialAccelerationZ = 4.0;    // m/s2

        double deltaT=1.0;      // sec
        
        double processErrorInX = 20;   //process errors in process covariance matrix not too big not too small   m
        double processErrorInY = 10;   //process errors in process covariance matrix not too big not too small   m
        double processErrorInZ = 30;   //process errors in process covariance matrix not too big not too small   m
     
        double processErrorInVx = 5;   //process errors in process covariance matrix not too big not too small   m/s
        double processErrorInVy = 5;   //process errors in process covariance matrix not too big not too small   m/s
        double processErrorInVz = 5;   //process errors in process covariance matrix not too big not too small   m/s

        double observationErrorInX = 25.0; // m
        double observationErrorInY = 25.0; // m
        double observationErrorInZ = 25.0; // m

        double observationErrorInVx = 6.0;  // m/s
        double observationErrorInVy = 6.0;  // m/s
        double observationErrorInVz = 6.0;  // m/s

        double ObservedX;
        double ObservedY;
        double ObservedZ;

        double ObservedVx;
        double ObservedVy;
        double ObservedVz;

        TrackerMd myTracker = new TrackerMd(initialPositionX,initialPositionY,initialPositionZ,
                                            intitialVelocityX,intitialVelocityY,intitialVelocityZ,
                                            initialAccelerationX,initialAccelerationY,initialAccelerationZ, 
                                            deltaT, 
                                            processErrorInX, processErrorInY, processErrorInZ,
                                            processErrorInVx, processErrorInVy,processErrorInVz,
                                            observationErrorInX,observationErrorInY,observationErrorInZ, 
                                            observationErrorInVx, observationErrorInVy,observationErrorInVz);

        try{
           File file =  new File ("observations"); 
           System.out.println(file.getAbsolutePath());
           Scanner scan = new Scanner(file);
           scan.useLocale(Locale.US);
           int counter = 1;
           while (scan.hasNextDouble()){

               ObservedX  = scan.nextDouble();
               ObservedY  = scan.nextDouble();
               ObservedZ  = scan.nextDouble();
               ObservedVx = scan.nextDouble();
               ObservedVy = scan.nextDouble();
               ObservedVz = scan.nextDouble();

               myTracker.trackingStep(ObservedX, ObservedY, ObservedZ, 
                                      ObservedVx, ObservedVy, ObservedVz);

               System.out.println("\n***************   " + counter + ". Step values  ************************************* :\n");
               System.out.println("predicted state is    :");
               myTracker.getPredictedStateMatrix().printMatrix();
               System.out.println("predicted process covariance matrix is    :");
               myTracker.getPredictedProcessCovarianceMatrix().printMatrix();
               System.out.println("Kalman gain is    :");
               myTracker.getKalmanGain().printMatrix();
               System.out.println("observation is        :");
               myTracker.getObservationMatrix().printMatrix();
               System.out.println("current state is      :");
               myTracker.getCurrentState().printMatrix();
               System.out.println("updated process covariance matrix is      :");
               myTracker.getUpdatedProcessCovarianceMatrix().printMatrix();
               counter++;

           } 
            scan.close();       
        }catch (FileNotFoundException e1){
        e1.printStackTrace();
        }
    }
}
