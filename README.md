# gnss-lidar_calibrate
calibrateHandEye，AX=BX, interior point way

THIS SCRIPT IS USING IN CALI THE GNSS/LIDAR INSTATION PARAMENTERS.

AUTHER: xiecong

DATE:2022.12.20

DEPENDENCES:

CVXPY(CONVEX SOLVER)-CVX-A64

RUN WAY:

    1. Running slerp.mlx: Getting the relative pose of gnss and lidarodometry, slerp the gnss pose as well.
    
    2. Running Data_prepocess.mlx: Processing the data for correcting the error. Getting the A and B of AX=BX.
    
    3. There are two way of sloving AX=BX:
    
        a. Using interior point way: Running interior_point_cali.mlx
        
        b. Using IRM to sloving the QCQP(transforming the AX=XB to QCQP model): Running qcqp_cali.m
