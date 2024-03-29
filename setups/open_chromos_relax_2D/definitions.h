#define  PHYSICS                        MHD
#define  DIMENSIONS                     2
#define  COMPONENTS                     2
#define  GEOMETRY                       CARTESIAN
#define  BODY_FORCE                     VECTOR
#define  FORCED_TURB                    NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  CHARACTERISTIC_TRACING
#define  DIMENSIONAL_SPLITTING          NO
#define  NTRACER                        0
#define  USER_DEF_PARAMETERS            16

/* -- physics dependent declarations -- */

#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  DIVB_CONTROL                   DIV_CLEANING
#define  BACKGROUND_FIELD               NO
#define  AMBIPOLAR_DIFFUSION            NO
#define  RESISTIVITY                    NO
#define  HALL_MHD                       NO
#define  THERMAL_CONDUCTION             SUPER_TIME_STEPPING
#define  VISCOSITY                      NO
#define  ROTATING_FRAME                 NO

/* -- user-defined parameters (labels) -- */

#define  HALF_LOOP_L                    0
#define  DELTA_CH                       1
#define  T_CH                           2
#define  T0_APEX                        3
#define  NE_CH                          4
#define  TC_L09                         5
#define  STEP_R                         6
#define  STEP_b                         7
#define  STEP_rT                        8
#define  STEP_rprs                      9
#define  Bz0                            10
#define  VRW_DT1                        11
#define  VRW_DT2                        12
#define  VRW_AV_FULLDOM_MIN             13
#define  VRW_AV_LAYER_MIN               14
#define  VRW_X_MAX_LAYER                15

/* [Beg] user-defined constants (do not change this line) */

#define  UNIT_LENGTH                    1e8
#define  UNIT_VELOCITY                  sqrt(CONST_kB/(0.5*CONST_mp)*1e6)
#define  UNIT_DENSITY                   1e-15
#define  VTK_TIME_INFO                  YES
#define  GLM_EXTENDED                   YES
#define  INITIAL_SMOOTHING              YES
#define  CHAR_LIMITING                  YES
#define  LIMITER                        MC_LIM
#define  PRINT_TO_FILE                  NO
#define  H_MASS_FRAC                    1.0
#define  He_MASS_FRAC                   0.0
#define  INTERNAL_BOUNDARY              YES

/* user-defined code configuration labels (options) */
#define  CH_OUTFLOW                     0
#define  CH_STRATIFICATION              1
#define  CH_EXTRAPOLATION               2
#define  INIT_DOMAIN_SOLVE_HS_EQUIL     3
#define  INIT_DOMAIN_LOAD_HS_EQUIL      4

/* user-defined code configuration labels (definitions) */
#define  CH_BOUNDARY_PRS_RHO            CH_STRATIFICATION
#define  CH_BOUNDARY_B                  CH_EXTRAPOLATION

/* [End] user-defined constants (do not change this line) */
