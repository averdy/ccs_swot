#include "OBSFIT_OPTIONS.h"
C--  File obsfit_read_obs.F: Routines to read OBSFIT input files

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: OBSFIT_READ_OBS

C !INTERFACE:
      SUBROUTINE OBSFIT_READ_OBS(
     O                      fNb, 
     I                      vNb, 
     I                      irec,
     I                      vec_loc,
     I                      myThid )

C !DESCRIPTION:
C Reads an observation and its weight from a netcdf obsfit 
C input file

C !USES:
      IMPLICIT NONE
C == Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#ifdef ALLOW_OBSFIT
# include "netcdf.inc"
# include "OBSFIT_SIZE.h"
# include "OBSFIT.h"
#endif

C !INPUT/OUTPUT PARAMETERS:
C fNb, vNb  :: 
C irec      :: 
C myThid    :: my Thread Id number
C vec_loc   :: 
      INTEGER fNb, vNb
      INTEGER irec, myThid
      _RL     vec_loc
CEOP

#ifdef ALLOW_OBSFIT

C !LOCAL VARIABLES:
c tmpObsNo :: number of obs in the file
      INTEGER k, kk, err, varid1, tmpObsNo
      INTEGER vec_start,vec_count
      _RL     vec_tmp1(1000),vec_tmp2(1000)
      CHARACTER*(max_len_mbuf) msgbuf

C     ==================================================================

      IF ( (irec.LT.obsfit_minind_buff).OR.
     &  (irec.GT.obsfit_maxind_buff).OR.
     & (obsfit_curfile_buff.NE.fNb) ) THEN
       err = NF_INQ_DIMID(fiddata_obs(fNb),'iOBS', varid1)
       err = NF_INQ_DIMLEN(fiddata_obs(fNb), varid1, tmpObsNo)

       IF (obsfit_curfile_buff.NE.fNb) THEN
c no asumption on whether a forward or a backward loop is calling
        obsfit_minind_buff=max(1,irec-500+1)
        obsfit_maxind_buff=min(tmpObsNo,irec+500)
       ELSEIF (irec.LT.obsfit_minind_buff) THEN
c implies that a backward loop is calling
        obsfit_minind_buff=max(1,irec-999)
        obsfit_maxind_buff=irec
       ELSE
c implies that a forward loop is calling
        obsfit_minind_buff=irec
        obsfit_maxind_buff=min(tmpObsNo,irec+999)
       ENDIF

#ifdef ALLOW_DEBUG
       IF ( debugLevel .GE. debLevD ) THEN
        WRITE(msgbuf,'(a,5I9)')
     &    'buffer readobsfile ',
     &    obsfit_minind_buff, obsfit_maxind_buff,
     &    irec, ObsNo(fNb), tmpObsNo
        CALL print_message(
     &    msgbuf, standardmessageunit, SQUEEZE_RIGHT, myThid)
       ENDIF
#endif

       vec_start=obsfit_minind_buff
       vec_count=obsfit_maxind_buff-obsfit_minind_buff+1

         err = NF_INQ_VARID(fiddata_obs(fNb),obsfit_names,
     &    varid1 )
         err = NF_GET_VARA_DOUBLE(fiddata_obs(fNb), varid1, vec_start,
     &    vec_count, vec_tmp1)
         err = NF_INQ_VARID(fiddata_obs(fNb),obsfit_namesweight
     &    , varid1 )
         err = NF_GET_VARA_DOUBLE(fiddata_obs(fNb), varid1, vec_start,
     &    vec_count, vec_tmp2)

       IF (err.NE.NF_NOERR) THEN
          WRITE(errorMessageUnit,'(A)')
     &     'WARNING in obsfit_read_obs: record not found!!'
       ENDIF

        DO kk=1,vec_count
          obsfit_data_buff(kk)=vec_tmp1(kk)
          obsfit_weight_buff(kk)=vec_tmp2(kk)
        ENDDO

       obsfit_curfile_buff=fNb
      ENDIF

cgf ...now, get vec_loc from the buffer
      IF (vNb.LT.0) THEN
       vec_loc=obsfit_weight_buff
     &   (irec-obsfit_minind_buff+1)

      ELSE
       vec_loc=obsfit_data_buff
     &   (irec-obsfit_minind_buff+1)
      ENDIF

#endif

      RETURN
      END
