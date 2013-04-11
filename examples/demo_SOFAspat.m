SOFAstart;
Obj=SOFAload('HRTF ARI NH2.sofa');
in=randn(5*Obj.Data.SamplingRate,1);
[out,azi,ele]=SOFAspat(in,Obj,[-45 45 0 15],[0 -30 90]);
subplot(2,1,1); plot(azi);
subplot(2,1,2); plot(ele);
wavplay(out,Obj.Data.SamplingRate);