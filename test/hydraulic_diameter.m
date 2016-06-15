%Calculate the diameter of the channels

%Given diameter of the chanell in mm
height = 1.6;   
length = 10;

%Calculate the radius of the arc
r = ((length/2)^2+height^2)/(2*height);

%Calculate the area
theta = asin(length/(2*r));
area = (theta*r^2)-0.5*length*(r-height);

%Calculate the perimeter
perimeter = 2*theta*r+length;

%Calculate the hydraulic diameter
D_H = 4*area/perimeter

a = 0.5*(perimeter/2+sqrt((perimeter/2)^2-4*area));
b = 0.5*(perimeter/2-sqrt((perimeter/2)^2-4*area));