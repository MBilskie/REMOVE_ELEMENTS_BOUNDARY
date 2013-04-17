!   MATTHEW BILSKIE
!   01/02/2013
!   DETERMINE IF A POINT IS OUTSIDE/INSIDE A POLYGON
!   PROGRAM MODIFIED FROM (http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html#Fortran Code for the Point in Polygon Test)
!
!   V3 - READ A POLYGON SHAPEFILE INSTEAD OF ASCII TEXT FILE
!
!   ----------------------------------------------------------------------------
!   VARIABLE DECLARATIONS
!   PX  -   X-COORDIANTE OF POINT IN QUESTION.
!   PY  -   Y-COORDINATE OF POINT IN QUESTION.
!   XX  -   N LONG VECTOR CONTAINING X-COORDINATES OF VERTICIES OF POLYGON.
!   YY  -   N LONG VECTOR CONTAINING Y-COORDINATES OF VERTICIES OF POLYGON.
!   N   -   NUMBER OF VERTICIES IN THE POLYGON.
!   INOUT   -   THE SIGNAL RETURNED:
!               -1  IF POINT IS OUTSIDE OF THE POLYGON.
!                0  IF POINT IS ON AN EDGE OR AT A VERTEX.
!                1  IF POINT IS INSIDE OF THE POLYGON.
!   ----------------------------------------------------------------------------
!   THE BOUNDARY INPUT FILE CONTAINS THE NUMBER OF VERTICIES AND THEN THE X Y
!   LOCATION OF EACH VERTEX.  THE VERTICIES CAN BE EITHER IN CCW OR CW ORDER.
!
!   THIS CAN BE CREATED IN SMS BY CLICKING OUT A CLOSED POLYLINE.  CONVERT THE
!   MAP COVERAGE TO A SCATTER SET (WITHOUT TRIANGULATING) AND EXPORT TO AN
!   ASCII TEXT FILE.
!   ----------------------------------------------------------------------------

    program PNPOLY
    
        character*60 in_mesh,out_mesh,in_boundary
        integer N,i,j,k
        integer NE,NP,NHY
        integer counter,newElem,newNodes
        real*8 PX,PY
        real*8 t1,t2
        logical MX,MY,NX,NY
        integer, allocatable :: INOUT(:)
        integer, allocatable :: NID(:),EID(:),NM(:,:),removeNodes(:)
        real*8, allocatable :: XX(:),YY(:),X(:),Y(:)
        real*8, allocatable :: xm(:),ym(:),zm(:)
        
        write(*,*)'NAME OF BOUNDARY FILE?'
        READ(*,*)in_boundary
        write(*,*)
        write(*,*)'NAME OF INPUT MESH'
        read(*,*)in_mesh
        write(*,*)
        write(*,*)'NAME OF OUTPUT MESH'
        read(*,*)out_mesh
        write(*,*)
        open(unit=100,file=in_boundary,status='old')
        read(100,*)N
        allocate (XX(N))
        allocate (YY(N))
        allocate (X(N))
        allocate (Y(N))
        do i=1,N
            read(100,*)XX(i),YY(i)
        enddo
        close(100)
        
        open(unit=14,file=in_mesh,status='old')
        read(14,*)
        read(14,*)NE,NP
        
        write(*,*)N,' BOUNDARY POINTS ON POLYGON'
        write(*,*)NE,' ELEMENTS'
        write(*,*)NP, ' NODES'
        write(*,*)
        
        allocate (INOUT(NP))
        allocate (NID(NP))
        allocate (xm(NP))
        allocate (ym(NP))
        allocate (zm(NP))
        allocate (removeNodes(NP))
        allocate (EID(NE))
        allocate (NM(NE,3))
        do i=1,NP
            read(14,*)NID(i),xm(i),ym(i),zm(i)
            if (i.ne.(NID(i))) STOP
        enddo
        do i=1,NE
            read(14,*)EID(i),NHY,NM(i,1),NM(i,2),NM(i,3)
        enddo
        close(14)
        write(*,*)'NODE/ELEMENT TABLE READ'
    
        counter=0
        write(*,*)'SEARCHING FOR NODES INSIDE/OUTSIDE POLYGON'
        write(*,*)
        call cpu_time(t1)
        do k=1,NP
            PX=xm(k)
            PY=ym(k)
            
            do 1 i=1,N
                X(i)=XX(i)-PX
1               Y(i)=YY(i)-PY
            INOUT(k)=-1
            do 2 i=1,N
                j=1+MOD(i,N)
                MX=X(i).GE.0.0
                NX=X(j).GE.0.0
                MY=Y(i).GE.0.0
                NY=Y(j).GE.0.0
                if (.not.((MY.or.NY).and.(MX.OR.NX)).or.(MX.and.NX)) GOTO 2
                if (.not.(MY.and.NY.and.(MX.or.NX).and..not.(MX.and.NX))) GOTO 3
                INOUT(k)=-INOUT(k)
                GOTO 2
3               IF ((Y(i)*X(j)-x(i)*y(j))/(X(j)-X(i))) 2,4,5
4               INOUT(k)=0
                GOTO 25
5               INOUT(k)=-INOUT(k)
2               continue
    
            if (INOUT(k).lt.0.0) then ! POINT IS OUTSIDE THE POLYGON
                counter=counter+1
                removeNodes(counter)=NID(k)
            endif
25      continue
        enddo ! end k loop
!25      continue

        write(*,*)'REMOVING OUTSIDE ELEMENTS...'
        write(*,*)
        newNodes=NP-counter
        newElem=0
        do i=1,NE
            !write(*,*)i,NE
            do j=1,counter
                if ( (NM(i,1).eq.removeNodes(j)).or.(NM(i,2).eq.removeNodes(j)).or.(NM(i,3).eq.removeNodes(j))) then
                    GOTO 35
                endif
            enddo
            newElem=newElem+1
35          continue
        enddo

        write(*,*)'WRITING NEW MESH FILE...'
        open(144,file=out_mesh,status='unknown')
        write(144,*)'REMOVED ELEMENTS'
        write(144,*)newElem,newNodes
        do i=1,NP
            if ( (INOUT(i).eq.0.0).or.(INOUT(i).gt.0.0) ) then
                write(144,*)NID(i),xm(i),ym(i),zm(i)
            endif 
        enddo

        do i=1,NE
            do j=1,counter
                if ( (NM(i,1).eq.removeNodes(j)).or.(NM(i,2).eq.removeNodes(j)).or.(NM(i,3).eq.removeNodes(j))) then
                    GOTO 45
                endif
            enddo
            write(144,*)EID(i),NHY,NM(i,1),NM(i,2),NM(i,3)
45          continue
        enddo
        
        write(*,*)
        call cpu_time(t2)
        !write(*,*)'Time in minutes =',(t2-t1)/60.
        
        close(144)
      
    end program