matlab_data = load('matlab_save_data.mat','data');
matlab = matlab_data.data;

data_size = size(matlab, 2);
if data_size == 4
    num_body = 1;
    adams = load('body1.txt');
elseif data_size == 7
    num_body = 2;
    adams = load('body2.txt');
elseif data_size == 10
    num_body = 3;
    adams = load('body3.txt');
elseif data_size == 13
    num_body = 4;
    adams = load('body4.txt');
elseif data_size == 16
    num_body = 5;
    adams = load('body5.txt');
elseif data_size == 19
    num_body = 6;
    adams = load('body6.txt');
end
index = 2;

ylabel_text = {'Position [rad]', 'Velocity [rad/s]', 'Acceleration [rad/s^2]'};

for i = 1 : num_body
    figure
    set(gcf,'color',[1,1,1])
    for j = 1 : 3
        subplot(3, 1, j)
        plot(matlab(:,1), matlab(:,index), 'LineWidth', 2)
        hold on
        plot(adams(:,1), adams(:,index),'--','LineWidth',2)
        grid on
        xlabel('Time [s]')
        ylabel(ylabel_text(j))
        set(gca,'FontSize',13)
        legend('MATLAB','ADAMS')
        index = index + 1;
    end
end