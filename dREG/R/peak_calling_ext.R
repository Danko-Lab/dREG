# First derivative of vector using 2-point central difference. T. C. O'Haver, 1988.
deriv<-function(a)
{
  n = NROW( a );
  d = rep(0, n );
  d[1] = a[2]-a[1];
  d[n] = a[n]-a[n-1];
  for (j in 2:(n-1))
      d[j] = ( a [j+1] - a[j-1] )/2;

  return(d);
}

#   SegmentedSmooth(y,w,type,ends) divides y into a number of equal-length
#   segments defined by the length of the vector 'smoothwidths', then
#   smooths each segment with smooth of type 'type' and width defined by
#   the elements of vector 'smoothwidths'.
#    Version 1, November 2016.
#   The argument "type" determines the smooth type:
#     If type=1, rectangular (sliding-average or boxcar)
#     If type=2, triangular (2 passes of sliding-average)
#     If type=3, pseudo-Gaussian (3 passes of sliding-average)
#     If type=4, pseudo-Gaussian (4 passes of same sliding-average)
#     If type=5, multiple-width (4 passes of different sliding-averages)
#   The argument "ends" controls how the "ends" of the signal
#   (the first w/2 points and the last w/2 points) are handled.
#     If ends=0, the ends are zero. (In this mode the elapsed
#       time is the fastest and is independent of the smooth width).
#     If ends=1, the ends are smoothed with progressively
#       smaller smooths the closer to the end. (In this mode the
#       elapsed time increases with increasing smooth widths).
#   SegmentedSmooth(Y,w,type) smooths with ends=0.
#   SegmentedSmooth(Y,w) smooths with type=1 and ends=0.
#
#   Examples: 3-segment smooth of random white noise, smooth widths of
#   2,20, and 200.
#     x=1:10000;y=randn(size(x));
#     plot(x,SegmentedSmooth(y,[2 20 200],3,0))
#
#     20-segment smooth, odd smooth widths from 1 to 41:
#     plot(x,SegmentedSmooth(y,[1:2:41],3,0))
#
#   Copyright (c) 2012, Thomas C. O'Haver

SegmentedSmooth<-function(y, smoothwidths, type=1, ends=0)
{
  NumSegments <- NROW(smoothwidths);
  SegLength <- round( NROW(y)/NumSegments );
  SmoothSegment <- matrix(0, nrow=NROW(y),ncol=NumSegments);
  SmoothedSignal <- rep(0, ncol=NROW(y));
  for (Segment in 1:NumSegments)
  {
    SmoothSegment[,Segment] <- fastsmooth( y, smoothwidths[Segment], type, ends);
      startindex = ( 1 + ( Segment-1)*SegLength);
      endindix = startindex + SegLength - 1;

      if (endindix > NROW(y) ){
      endindix <- NROW(y)
        break;
    }

      SmoothedSignal[startindex:endindix] <- SmoothSegment[startindex:endindix, Segment];
  }

  SmoothedSignal;
}


# fastsmooth(Y,w,type,ends) smooths vector Y with smooth
#  of width w. Version 3.0, October 2016.
# The argument "type" determines the smooth type:
#   If type=1, rectangular (sliding-average or boxcar)
#   If type=2, triangular (2 passes of sliding-average)
#   If type=3, pseudo-Gaussian (3 passes of sliding-average)
#   If type=4, pseudo-Gaussian (4 passes of same sliding-average)
#   If type=5, multiple-width (4 passes of different sliding-average)
# The argument "ends" controls how the "ends" of the signal
# (the first w/2 points and the last w/2 points) are handled.
#   If ends=0, the ends are zero.  (In this mode the elapsed
#     time is independent of the smooth width). The fastest.
#   If ends=1, the ends are smoothed with progressively
#     smaller smooths the closer to the end. (In this mode the
#     elapsed time increases with increasing smooth widths).
# fastsmooth(Y,w,type) smooths with ends=0.
# fastsmooth(Y,w) smooths with type=1 and ends=0.
# Examples:
# fastsmooth([1 1 1 10 10 10 1 1 1 1],3)= [0 1 4 7 10 7 4 1 1 0]
#
# fastsmooth([1 1 1 10 10 10 1 1 1 1],3,1,1)= [1 1 4 7 10 7 4 1 1 1]
#
# x=1:100;
# y=randn(size(x));
# plot(x,y,x,fastsmooth(y,5,3,1),'r')
# xlabel('Blue: white noise.    Red: smoothed white noise.')
#
# Copyright (c) 2012, Thomas C. O'Haver
#

fastsmooth<-function( Y, w, type, ends )
{
    if( type == 1 )
       SmoothY=sa(Y,w,ends)
    else if( type == 2 )
       SmoothY=sa(sa(Y,w,ends),w,ends)
    else if( type == 3 )
       SmoothY=sa(sa(sa(Y,w,ends),w,ends),w,ends)
    else if( type == 4 )
       SmoothY=sa(sa(sa(sa(Y,w,ends),w,ends),w,ends),w,ends)
    else if( type == 5 )
       SmoothY=sa(sa(sa(sa(Y,round(1.6*w),ends),round(1.4*w),ends),round(1.2*w),ends),w,ends)

  return(SmoothY);
}

sa<-function(Y, smoothwidth, ends)
{
  w <- round(smoothwidth);

  SumPoints <- sum(Y[1:w]);
  s <- rep( 0, NROW(Y) );
  halfw <- round(w/2);

  L <- NROW(Y);
  for (k in 1:(L-w))
  {
    s[k+halfw-1] <- SumPoints;
      SumPoints <- SumPoints - Y[k];
      SumPoints <- SumPoints + Y[k+w];
  }

  s[ k+halfw] = sum ( Y[(L-w+1):L] );
  SmoothY <- s/w;

  # Taper the ends of the signal if ends=1.
   if (ends==1)
    {
    startpoint <- (smoothwidth + 1)/2;
      SmoothY[1]=( Y[1] + Y[2])/2;
      for(k in 2:startpoint )
      {
      SmoothY[k] = mean(Y [1:(2*k-1)]);
         SmoothY[L-k+1] = mean(Y[(L-2*k+2):L]);
    }

      SmoothY[L]=( Y[L] + Y[L-1] )/2;
  }

  SmoothY
}
