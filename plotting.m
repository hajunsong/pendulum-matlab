% matlab_data = load('matlab_save_data.mat','data');
% matlab = matlab_data.data;

% data_size = size(matlab, 2);
% if data_size == 4
%     num_body = 1;
%     adams = load('body1.txt');
% elseif data_size == 7
%     num_body = 2;
%     adams = load('body2.txt');
% elseif data_size == 10
%     num_body = 3;
%     adams = load('body3.txt');
% elseif data_size == 13
%     num_body = 4;
%     adams = load('body4.txt');
% elseif data_size == 16
%     num_body = 5;
%     adams = load('body5.txt');
% elseif data_size == 19
%     num_body = 6;
%     adams = load('body6.txt');
% end
num_body = 1;
matlab = load(sprintf('matlab_body%d.txt',num_body));
adams = load(sprintf('body%d.txt',num_body));
index = 2;

ylabel_text = {'Position [rad]', 'Velocity [rad/s]', 'Residual [Nm]'};

for i = 1 : num_body
    figure
    set(gcf,'color',[1,1,1])
    for j = 1 : 3
        subplot(3, 1, j)
%         if index == 4
%             plot(matlab(:,1), matlab(:,index+1), 'LineWidth', 2)
%         else
            plot(matlab(:,1), matlab(:,index), 'LineWidth', 2)
%         end
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

figure
set(gcf,'color',[1,1,1])
subplot(311)
plot(matlab(:,1), matlab(:,11), 'LineWidth',2)
grid on
xlabel('Time [s]')
ylabel('Rasidual [Nm]')

subplot(312)
plot(matlab(:,1), matlab(:,12), 'LineWidth',2)
grid on
xlabel('Time [s]')
ylabel('Rasidual [Nm]')

subplot(313)
plot(matlab(:,1), matlab(:,13), 'LineWidth',2)
grid on
xlabel('Time [s]')
ylabel('Rasidual [Nm]')
set(gca,'FontSize',13)