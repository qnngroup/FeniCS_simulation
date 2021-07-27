%template with the right precision setup
template = 'C:\Users\John\Desktop\Lorentz_Scripted\solution.pvd';

%variables initialization
iErr = int16(1);

Lorentz = actxserver('IES.document');
iErr = 1;
while iErr > 0      % pause until Lorentz is completely initialized
  [iErr] = Window_CheckInitialization(Lorentz, iErr);  
end
[iErr] = File_OpenTemplate(Lorentz,template, iErr);     %load template
%Model_Delete_All(Lorentz); %clear everything to begin