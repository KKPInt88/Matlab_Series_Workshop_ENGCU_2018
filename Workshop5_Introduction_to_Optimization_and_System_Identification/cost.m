function J = cost(X, XZ_ref)

    % Given L12;
    Ts = 0.001;
    L1 = X(1);
    L2 = X(2);
    
    
%     path_Connector = 'Part2_5_Four_Bar2/Connector Link';
%     path_Rocker    = 'Part2_5_Four_Bar2/Rocker Link';  %gcb(mdl)
%     
%     % Assign values in base workspace
%     assignin('base','a',param_v(1));
%     assignin('base','b',param_v(2));
%     
%     % Assign values in base workspace
%     assignin('Connector Link','L',L1);
%     assignin('Rocker Link','L',L2);

    %get_param(gcb,'MaskValueString')
    %set_param(path_Connector, 'MaskValueString', sprintf('1000|%.0f|2|0.5|0.5|[0.6 0.6 0.0]|sm_four_bar_connector.png',L1) ); 
    %set_param(path_Rocker,    'MaskValueString', sprintf('1200|%.0f|2|0.5|0.5|[1.0 0.2 0.2]|sm_four_bar_rocker.png',L2)    );
    
    simopt = simset('SrcWorkspace','Current');
    set_param(gcs,'SimulationCommand','Update'); % Update Model
    
    [~, ~, Yout]= sim('Part2_6_Four_Bar', [0:1:3000]'*Ts, simopt);
    
    x = Yout(:,1);
    z = Yout(:,2);

    x_ref = XZ_ref(:,1);
    z_ref = XZ_ref(:,2);

    J = sum( (x - x_ref).^2 + (z - z_ref).^2 );
    
    % Plot 
    figure(101);
%     plot(x_ref, z_ref, 'LineWidth', 2, 'Color', 'b');
%     hold on;
    plot(x, z, 'LineWidth', 4, 'LineStyle', '--');
    
 
end