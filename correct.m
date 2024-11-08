clear all;

% Načtení dat z Excelu
data = readtable('compiled.xlsx');
sphere = readmatrix('sphere_wall_reflectance3.csv');

% Vynechání prvního sloupce (vlnové délky)
wavelengths = data{2:end, 1}; % První sloupec, třetí řádek a dále
gapFraction = data{1, 2:end}; % Hodnoty gap fraction od druhého sloupce ve druhém řádku
gapFraction2 = gapFraction(:, 1:2:end);

% Inicializace proměnných pro opravené reflektance a transmittance
Reflectance = data{2:end, 2:2:end}; % R_ sloupce bez prvního sloupce
Transmittance = data{2:end, 3:2:end}; % T_ sloupce bez prvního sloupce

[~, idx_799] = min(abs(wavelengths - 799)); % Najdeme nejbližší hodnotu 799 v wavelengths
[~, idx_550] = min(abs(wavelengths - 550)); % Najdeme nejbližší hodnotu 550 v wavelengths
[~, idx_450] = min(abs(wavelengths - 450)); % Najdeme nejbližší hodnotu 450 v wavelengths

alternativeGF_R1 = (0.49 - Reflectance(idx_799, :)) / 0.58; % Výpočet alternativních GF
alternativeGF_R2 = (0.16 - Reflectance(idx_550, :)) / 0.19; % Výpočet alternativních GF
alternativeGF_R3 = (0.06 - Reflectance(idx_450, :)) / 0.06; % Výpočet alternativních GF

alternativeGF_T1 = (Transmittance(idx_799, :) - 0.34) / 0.59; % Výpočet alternativních GF
alternativeGF_T2 = (Transmittance(idx_550, :) - 0.09) / 0.88; % Výpočet alternativních GF
alternativeGF_T3 = (Transmittance(idx_450, :) - 0.03) / 0.93; % Výpočet alternativních GF

% Mesarch et al. korekce pro reflektanci a transmittanci
for i = 1:size(Reflectance, 2)
    correctedReflectance(:, i) = Reflectance(:, i) ./ (1 - alternativeGF_R1(i));
    correctedReflectance2(:, i) = Reflectance(:, i) ./ (1 - alternativeGF_R2(i));
    correctedReflectance3(:, i) = Reflectance(:, i) ./ (1 - alternativeGF_R3(i));
    correctedReflectance4(:, i) = Reflectance(:, i) ./ (1 - gapFraction2(i));
    
    correctedTransmittance(:, i) = (Transmittance(:, i) - (sphere(:,2)*alternativeGF_T1(i))) ./ (1 - alternativeGF_T1(i));
    correctedTransmittance2(:, i) = (Transmittance(:, i) - (sphere(:,2)*alternativeGF_T2(i))) ./ (1 - alternativeGF_T2(i));
    correctedTransmittance3(:, i) = (Transmittance(:, i) - (sphere(:,2)*alternativeGF_T3(i))) ./ (1 - alternativeGF_T3(i));
    correctedTransmittance4(:, i) = (Transmittance(:, i) - (sphere(:,2)*alternativeGF_R1(i))) ./ (1 - alternativeGF_R1(i));
    correctedTransmittance5(:, i) = (Transmittance(:, i) - (sphere(:,2)*gapFraction2(i))) ./ (1 - gapFraction2(i));
end

% Define colors for Reflectance and 1 - Transmittance for each sample
reflectanceColor = [0.2, 0.4, 0.6]; % Blueish color for Reflectance
transmittanceColor = [0.8, 0.2, 0.2]; % Reddish color for 1 - Transmittance

% Indices for every 10th wavelength to display error bars
errorBarIndices = 1:10:length(wavelengths);

% Plot average corrected reflectance and 1 - average corrected transmittance for each sample set with selected error bars
figure;
for i = 1:6
    subplot(3, 2, i); % Arrange plots in a 3x2 grid

    % Select the corresponding corrected reflectance and transmittance data, calculate averages and standard deviations
    switch i
        case 1
            avgReflectance = mean(correctedReflectance, 2);
            stdReflectance = std(correctedReflectance, 0, 2);
            avgTransmittance = mean(correctedTransmittance, 2);
            stdTransmittance = std(correctedTransmittance, 0, 2);
            titleText = 'Sample Set 1';
        case 2
            avgReflectance = mean(correctedReflectance2, 2);
            stdReflectance = std(correctedReflectance2, 0, 2);
            avgTransmittance = mean(correctedTransmittance2, 2);
            stdTransmittance = std(correctedTransmittance2, 0, 2);
            titleText = 'Sample Set 2';
        case 3
            avgReflectance = mean(correctedReflectance3, 2);
            stdReflectance = std(correctedReflectance3, 0, 2);
            avgTransmittance = mean(correctedTransmittance3, 2);
            stdTransmittance = std(correctedTransmittance3, 0, 2);
            titleText = 'Sample Set 3';
        case 4
            avgReflectance = mean(correctedReflectance4, 2);
            stdReflectance = std(correctedReflectance4, 0, 2);
            avgTransmittance = mean(correctedTransmittance4, 2);
            stdTransmittance = std(correctedTransmittance4, 0, 2);
            titleText = 'Sample Set 4';
        case 5
            avgReflectance = mean(correctedReflectance4, 2); % Assuming correctedReflectance4 is intended
            stdReflectance = std(correctedReflectance4, 0, 2);
            avgTransmittance = mean(correctedTransmittance5, 2);
            stdTransmittance = std(correctedTransmittance5, 0, 2);
            titleText = 'Sample Set 5';
        case 6
            avgReflectance = mean(correctedReflectance, 2); % Assuming correctedReflectance4 is intended
            stdReflectance = std(correctedReflectance, 0, 2);
            avgTransmittance = mean(correctedTransmittance2, 2);
            stdTransmittance = std(correctedTransmittance2, 0, 2);
            titleText = 'Sample Set 6';
    end

    % Plot Average Corrected Reflectance
    plot(wavelengths, avgReflectance, 'Color', reflectanceColor, 'DisplayName', 'Average Reflectance');
    hold on;
    
    % Plot 1 - Average Corrected Transmittance
    plot(wavelengths, 1 - avgTransmittance, 'Color', transmittanceColor, 'DisplayName', '1 - Average Transmittance');
    
    % Plot selected error bars
    errorbar(wavelengths(errorBarIndices), avgReflectance(errorBarIndices), stdReflectance(errorBarIndices), ...
        'o', 'Color', reflectanceColor, 'MarkerSize', 4, 'LineWidth', 1.2, 'DisplayName', 'Reflectance Error Bars');
    
    errorbar(wavelengths(errorBarIndices), 1 - avgTransmittance(errorBarIndices), stdTransmittance(errorBarIndices), ...
        'o', 'Color', transmittanceColor, 'MarkerSize', 4, 'LineWidth', 1.2, 'DisplayName', '1 - Transmittance Error Bars');
    
    % Add title, labels, and legends
    title(titleText);
    xlabel('Wavelength (nm)');
    ylabel('Value');
    xlim([350, 2500]); % Set x-axis limits
    ylim([0, 1]);
    legend;
    hold off;
end

% Adjust layout
sgtitle('Average Corrected Reflectance and 1 - Corrected Transmittance with Selected Standard Deviation Error Bars');
