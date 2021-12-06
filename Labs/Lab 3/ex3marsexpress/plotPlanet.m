
function [globe] = plotPlanet(planetID, position, handle, scaleFactor)

    % Preraring Figure Object
    if nargin<3
        HAXIS = gca;
    elseif ishandle(handle)==0
            msg = ['The figure handle is not valid'];
            eid = sprintf('TOOLBOX:%s:propertyError', mfilename);
            error(eid,'%s',msg)
    else
        try
            HAXIS=gca(handle);
        catch
            HAXIS=handle;  
        end
        hold on
    end
    %--------------------------------------------------------------------------
    if nargin<4
        scale=1;   
    end
    
    planetNames = {'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto', 'Sun'};
    


    Rplanet = astroConstants(3)*scaleFactor; % Planet radius w.r.t. sun [km]
    npanels = 360; % Number of globe panels around the equator [deg/panel] = [360/npanels]
    alpha = 0.9;  % alpha (i.e. transparency level) of the globe
    erad=Rplanet; % equatorial radius [km]
    prad=Rplanet; % polar radius [km]
    hold on;
    axis equal;
    axis vis3d;
    % Create a 3D meshgrid of the sphere points using the ellipsoid function
    [x,y,z] = ellipsoid(position(1), position(2), position(3), erad, erad, prad, npanels);
    globe = surf(HAXIS, x,y,z,'FaceColor','none','EdgeColor',0.5*[1 1 1], 'HandleVisibility','off');
    % RMK.: HandleVisibility=off removes unuseful legends for the plotted
    % globe
    cdata=imread(sprintf('%s.jpg',planetNames{planetID})); % Load Earth image for texture map
    % Set image as color data (cdata) property, and set face color to indicate
    % a texturemap, which Matlab expects to be in cdata.

    globe.FaceColor = 'texturemap';
    globe.CData = cdata;
    globe.EdgeColor = 'none';
    hfig = gcf; 


%      globe.AmbientStrength = 0.1;
%      globe.DiffuseStrength = 1;
%      globe.SpecularColorReflectance = .5;
%      globe.SpecularExponent = 20;
%      globe.SpecularStrength = 1;
%      globe.FaceLighting = 'gouraud';

     globe.FaceLighting = 'gouraud';
     globe.AmbientStrength = 0.5;
     if planetID ~= 10
        globe.FaceAlpha = 0.9;
     else
        globe.FaceAlpha = 1.0;
        globe.FaceLighting = 'none';
     end
    

end