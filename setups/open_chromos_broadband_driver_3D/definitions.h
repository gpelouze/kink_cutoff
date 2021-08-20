#define  PHYSICS                        MHD
#define  DIMENSIONS                     3
#define  COMPONENTS                     3
#define  GEOMETRY                       CARTESIAN
#define  BODY_FORCE                     POTENTIAL
#define  FORCED_TURB                    NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  CHARACTERISTIC_TRACING
#define  DIMENSIONAL_SPLITTING          NO
#define  NTRACER                        0
#define  USER_DEF_PARAMETERS            7

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
#define  TC_L09                         1
#define  STEP_R                         2
#define  STEP_b                         3
#define  DRIVER_v0                      4
#define  VRW_AV_LAYER_MIN               5
#define  VRW_X_MAX_LAYER                6

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
#define  DEBUG                          NO

/* user-defined code configuration labels (options) */
#define  CH_OUTFLOW                     0
#define  CH_STRATIFICATION              1
#define  CH_EXTRAPOLATION               2

/* user-defined code configuration labels (definitions) */
#define  CH_BOUNDARY_PRS_RHO            CH_STRATIFICATION
#define  CH_BOUNDARY_B                  CH_EXTRAPOLATION

/* [End] user-defined constants (do not change this line) */
