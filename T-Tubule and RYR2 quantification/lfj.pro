pro lfj
;***********for images recorded by confocal 510

;loadct,33

  close,2
  s=strarr(51)
  s(1)='a'&s(2)='b'&s(3)='c'&s(4)='d'&s(5)='e'&s(6)='f'&s(7)='g'&s(8)='h'&s(9)='i'&s(10)='j'
  s(11)='k'&s(12)='l'&s(13)='m'&s(14)='n'&s(15)='o'&s(16)='p'&s(17)='q'&s(18)='r'&s(19)='s'&s(20)='t'
  s(21)='u'&s(22)='v'&s(23)='w'&s(24)='x'&s(25)='y'&s(26)='z'&s(27)='aa'&s(28)='ab'&s(29)='ac'&s(30)='ad'
  s(31)='ae'&s(32)='af'&s(33)='ag'&s(34)='ah'&s(35)='ai'&s(36)='aj'&s(37)='ak'&s(38)='al'&s(39)='am'&s(40)='an'
  s(41)='ao'&s(42)='ap'&s(43)='aq'&s(44)='ar'&s(45)='as'&s(46)='at'&s(47)='au'&s(48)='av'&s(49)='aw'&s(50)='ax'
  
  ;**************************  Input parameter settings  ***********************************************
	;openu,2,'z:\zhoupeng\result\�½� �ı��ĵ�.txt', /APPEND
	ser1=1
	ser2=1
	readdir1='G:\CMYA5 quantification-Myh6 KO\CMYA5 quantification'	        ;directory of file to read
	;readdir2='11'        ;directory of file to read
	black=0                  					;backgroud sbstraction1
	;Ca0=100										;resing Ca concentration
  for j=ser1,ser2 do begin
  
  	fileser=s(j)			   ;the starting letter of file series
    for	frstfil=56,56 do begin        ;begaining file number
      ;file=readdir1+'\Image4_Copychannel.lsm';+strtrim(string(frstfil),2)+'.lsm'
      file = readdir1+'\20210408-dish03-AMCMS-Myh6 casAAV-anti-CMYA5 555-anti-SAA 647-f09_C003T001.tif'
      if file_test(file) eq 1 then begin ;��֤�Ƿ�����ļ���1=�У���ʼ����
        filetype = strmid(file, strlen(file)-3, 3)
        if filetype eq 'lsm' then begin
          m=lonarr(2);�������ݸ���Ϊ2
          n=intarr(2)
          close,1
          openr,1,file
          readu,1,m								    ;get offset to the first image directory
          POINT_LUN,1,m[1]+2
          while (n[0] ne 256) do readu,1,n,m & w=m[1]	;image width (x-direction)
          while (n[0] ne 257) do readu,1,n,m & t=m[1]	;image length (lines)
          while (n[0] ne 273) do readu,1,n,m & p=m[1]	;offset of image
          while ((n[0] ne -31124) and (n[0] ne 34412)) do readu,1,n,m
          n=dblarr(1)
          image=bytarr(w,t)
          POINT_LUN,1,m[1]+40
          readu,1,n				;offset of LSM info
          ww=n[0]*1000000								;pixel distance in x-direction (in micrometer)
          POINT_LUN,1,m[1]+132
          readu,1,m			;offset to time stamps
          POINT_LUN,1,m[0]+8
          readu,1,n			;absolute time of the first line
          tt=n[0]
          POINT_LUN,1,m[0]+t*8
          readu,1,n			;absolute time of the last line
          tt=(n[0]-tt)*1000/(t-1)						;time interval between lines (in milisecond)
          POINT_LUN,1,p
          readu,1,image				;reading the image
          close,1
          ;print,file,':',w,'  X',t,'  image with',ww,'um X',tt,'ms pixel size'
        endif else begin
          image = read_tiff(file)
          szimg = size(image, /dim)
          ww = 0.2
          w = szimg[0]
          t = szimg[1]
        endelse
        
        ;------------- image display  --------------------------
        image=image+0.                                ;�����Ϊ�������飬���������Ĳ���
        ima=smooth(image,3,edge=1)                    ;ƽ��ͼƬ��3Ϊƽ����С��λ��edge=1���ϵĵ�Ҳƽ��
        window,2,xs=w,ys=t,title=file+' WINDOW 2' 	;for original normalized image only
        ;window,0,xs=t,ys=w,title=file+' WINDOW 0'		;for original normalized image overlayed with 2SD sparks
        wshow,2
        wset,2
        tvscl,ima*3<255									;display original
        
        
        ;*********************** Spectrum **********************
        ;cursor,left,bottom,3,/device
        ;cursor,right,top,3,/device
        ;imb=ima[left:right,bottom:top]
        ;imb=imb-mean(imb)
        ;spectrum=fltarr(right-left+1)
        ;spectrum[*]=abs((fft(imb))[*,0])
        
        
        ;*********************** TTb Detection *****************
        
        ima=ima-min(ima)
        
        cursor,x11,y11,3,/device                   ;cursor�ǵ���궯����device�Ǹ��ͼƬ��С���ص����ֵ
        cursor,x22,y22,3,/device
        cell_thrshld=mean(ima(x11:x22,y11:y22))
        saturation=max(ima)
        ;his=histogram(ima)
        window,1,xs=640,ys=512
        tvscl,ima
        cursor,x1,y1,3,/device & cursor,x2,y2,3,/device
        ;plot,his
        ;xyouts,300,400,'Please point out the cell_thrshld.',/device,charsize=2
        ;cursor,x,y,3
        ;temp=max(his[0:x],black)
        ;temp=min(his[black:x],cell_thrshld)
        ;temp=min(his[x:*],saturation)
        ;cell_thrshld=cell_thrshld+black
        ;saturation=saturation+x
        imc=fix(ima)*0+1
        imc[where(ima le cell_thrshld)]=0
        imc[where(ima ge saturation)]=0
        peak1=max(histogram(ima(x1:x2,y1:y2),min=0),middle)
        middle=middle+(max(ima(x1:x2,y1:y2))-middle)/3
        amplif=250/middle
        ima=ima*amplif+10<255
        imd=(ima/smooth(ima+.0,2/ww,edge=1))
        imd=imd-min(imd)
        imd=(imd+.0)/max(imd)
        tvscl,imd
        wset,1
        wshow,1
        imc=imc+0.
        imd=imd+0.
        rad=0
        repeat begin
        	tv,rot(imc*imd*200<255,rad,/interp,missing=0)			;display original
        	cursor,x_temp,y_temp,3
        	if !mouse.button eq 1 then rad=rad+10
        	if !mouse.button eq 4 then rad=rad-10
        ENDREP UNTIL !mouse.button eq 2
        
        repeat begin
        	tv,rot(imc*imd*200<255,rad,/interp,missing=0)			;display original
        	cursor,x_temp,y_temp,3
        	if !mouse.button eq 1 then rad=rad+1
        	if !mouse.button eq 4 then rad=rad-1
        ENDREP UNTIL !mouse.button eq 2
        imc1=imc
        imd1=imd
        lng=max([w,t])
        print,lng
        imc=fltarr(lng,lng)
        imd = imc
        imc[0:w-1,0:t-1]=imc1
        imd[0:w-1,0:t-1]=imd1
        imc=shift(imc,(lng-w)/2,(lng-t)/2)
        imd=shift(imd,(lng-w)/2,(lng-t)/2)
        
        imc=rot(imc,rad,/interp,missing=0)
        imd=rot(imd,rad,/interp,missing=0)
        window,2,xs=lng,ys=lng
        tv,imc*imd*200<255
        cursor,left,bottom,3,/device
        cursor,right,top,3,/device
        top = bottom + right - left
        imf=imd[left:right,bottom:top]
        ime=imc[left:right,bottom:top]
        imf=imf-min(imf)
        his=histogram(imf,bin=0.01)
        for ttb_thrshld=1,25500 do if total(his[1:ttb_thrshld]) ge total(his)*0.7 then break
        ttb_thrshld=ttb_thrshld*0.01
        ime[where(imf le ttb_thrshld)]=0
        
        ;*********************** Analysis **********************
        
        xx=right-left+1
        yy=top-bottom+1
        ;*******************��ȡ����**************************
        ;ime = (xx lt yy) ? ime[*,yy/2-xx/2:yy/2+xx/2] : ime[xx/2-yy/2:yy/2+xx/2,*]
        xx = (size(ime))[1]
        yy = (size(ime))[2]
        
        wdth=fix(2/ww)/2*2+1
        point=intarr(xx,yy,2)
        line=intarr(99,99,2)
        line[(44-wdth/2):(44+wdth/2),44,0]=1
        line[44,(44-wdth/2):(44+wdth/2),1]=1
        
        for k=0,1 do for angle=-20,20 do begin
        	sgmnt=rot(line[*,*,k],angle)
        	kx=fix(total(not(not(rebin(sgmnt+0.,99,1)))))/2*2+1
        	ky=fix(total(not(not(rebin(sgmnt+0.,1,99)))))/2*2+1
        	kernel=sgmnt[min(where(rebin(sgmnt+0.,99,1))):min(where(rebin(sgmnt+0.,99,1)))+kx-1,min(where(rebin(sgmnt+0.,1,99))):min(where(rebin(sgmnt+0.,1,99)))+ky-1]
        		conv=fix(convol(ime,kernel,total(kernel)))
        	for sub1=0+kx/2,xx-kx/2-1 do for sub2=0+ky/2,yy-ky/2-1 do if conv[sub1,sub2] then point[sub1-kx/2:sub1+kx/2,sub2-ky/2:sub2+ky/2,k]=point[sub1-kx/2:sub1+kx/2,sub2-ky/2:sub2+ky/2,k] or kernel
        endfor
        window,3,xs=xx,ys=yy,xp=0,yp=0
        tvscl,point[*,*,0]+ime
        window,4,xs=xx,ys=yy,xp=right-left+10,yp=0
        tvscl,point[*,*,1]+ime
        ;imfft=fft(ime)
        ;ime=ime+0.0
        imfft=abs(fft(ime/total(ime)))
        window,5
        plot,imfft[0,*]
        plot,rebin(imfft,1,yy)
        cursor,left,temp,3
        cursor,right,temp,3
        fft_peak=max(imfft[0,left:right],max_frequency)
        max_frequency=max_frequency+left
        print,file,'  ',total(point[*,*,0])/total(ime),'  ',total(point[*,*,1])/total(ime),'	',total(point[*,*,0])/total(ime*0+1),'	',total(point[*,*,1])/total(ime*0+1),'	',fft_peak,'	',(top-bottom+1)*ww/max_frequency
        ;�print,imfft[0,*]
        ;'  ',total(ime ne 0)/(xx*yy)
        ;�print,rebin(imfft,1,yy)
        ;write_path = file_dirname(file) + '/' + file_basename(file) + 'result.csv'
        ;write_csv,write_path,[total(point[*,*,0])/total(ime),total(point[*,*,1])/total(ime),total(point[*,*,0])/total(ime*0+1),total(point[*,*,1])/total(ime*0+1),fft_peak,(top-bottom+1)*ww/max_frequency]
      endif
    endfor
  endfor
  close,2
end

