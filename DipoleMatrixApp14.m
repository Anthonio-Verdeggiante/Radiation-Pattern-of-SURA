classdef DipoleMatrixApp14 < matlab.apps.AppBase

    % Свойства приложения
    properties (Access = public)
        UIFigure          matlab.ui.Figure
        GridLayout        matlab.ui.container.GridLayout % Сетка для разделения интерфейса
        PowerSlider1      matlab.ui.control.Slider % Мощность секции 1
        PowerSlider2      matlab.ui.control.Slider % Мощность секции 2
        PowerSlider3      matlab.ui.control.Slider % Мощность секции 3
        FrequencySlider   matlab.ui.control.Slider % Единый слайдер частоты
        ThetaSlider       matlab.ui.control.Slider % Слайдер для угла theta
        PhaseNoiseSlider  matlab.ui.control.Slider % Слайдер для СКО фазы
        AmplitudeNoiseSlider matlab.ui.control.Slider % Слайдер для СКО амплитуды
        RowDropDown       matlab.ui.control.DropDown % Выпадающий список для выбора ряда Nx
        ColumnDropDown    matlab.ui.control.DropDown % Выпадающий список для выбора столбца Ny
        PowerLabel1       matlab.ui.control.Label
        PowerLabel2       matlab.ui.control.Label
        PowerLabel3       matlab.ui.control.Label
        FrequencyLabel    matlab.ui.control.Label
        ThetaLabel        matlab.ui.control.Label  % Метка для угла theta
        PhaseNoiseLabel   matlab.ui.control.Label  % Метка для СКО фазы
        AmplitudeNoiseLabel matlab.ui.control.Label % Метка для СКО амплитуды
        RowLabel          matlab.ui.control.Label  % Метка для выбора ряда Nx
        ColumnLabel       matlab.ui.control.Label  % Метка для выбора столбца Ny
        ModeSwitch        matlab.ui.control.Switch % Тумблер для выбора режима (X-мода или O-мода)
        UIAxes            matlab.ui.control.UIAxes % Для интенсивности матрицы
        Matrix            % Матрица комплексных чисел (12x12x2)
        DirectionButton   matlab.ui.control.UIControl % Кнопка для отладки ЭМП
        DNPButton         matlab.ui.control.UIControl % Кнопка для построения ДН
        TextArea          matlab.ui.control.TextArea % Область для сообщений
        % Добавляем оси для графиков направленности
        AxesPsi           matlab.ui.control.UIAxes % График F_n(psi)
        AxesThetaX        matlab.ui.control.UIAxes % График F_n(theta_x)
        AxesThetaY        matlab.ui.control.UIAxes % График F_n(theta_y)
        AxesTheta3D       matlab.ui.control.UIAxes % График F_n(theta_x, theta_y)
        AxesXY            matlab.ui.control.UIAxes % График F_n(x, y)
        SelectedRow       double % Выбранный ряд Nx
        SelectedColumn    double % Выбранный столбец Ny
    end

    methods (Access = private)

        % Инициализация приложения
        function startupFcn(app)
            app.appendMessage('Инициализация приложения...');
            % Изначальная матрица 12x12
            app.SelectedRow = 1; % Начальный выбранный ряд
            app.SelectedColumn = 1; % Начальный выбранный столбец
            app.Matrix = initializeMatrix(app, [300, 300, 300], app.FrequencySlider.Value, app.ThetaSlider.Value, app.PhaseNoiseSlider.Value, app.AmplitudeNoiseSlider.Value);
            updatePlot(app); % Обновляем график
            app.appendMessage('График должен быть инициализирован');
            app.appendMessage(['']);
            % Первоначальное построение графиков направленности
            updateDirectionPlots(app);
        end

        % Функция добавления сообщения в TextArea
        function appendMessage(app, message)
            if isempty(app.TextArea.Value)
                app.TextArea.Value = {message};
            else
                app.TextArea.Value = [app.TextArea.Value; {message}];
            end
        end

        % Функция создания матрицы
        function matrix = initializeMatrix(app, powers, frequency, theta, phase_noise, amplitude_noise)
            matrix = zeros(12, 12, 2); % Матрица 12x12x2 для двух диполей на элемент
            c = 3e8; % Скорость света, м/с
            dipole_length = 33.830; % Длина диполя, м
            resonant_frequency = c / (2 * dipole_length); % Резонансная частота
            resonant_frequency_MHz = resonant_frequency / 1e6; % Перевод в МГц (~4.41 МГц)
            
            % Резонансная характеристика (гауссова функция)
            bandwidth = 5.0; % Полоса пропускания, МГц
            amplitude_factor = exp(-((frequency - resonant_frequency_MHz) / bandwidth).^2);
            
            % Шаг между диполями в метрах
            dx = 0.5; % Шаг в долях длины волны
            wavelength = c / (frequency * 1e6);
            dx_meters = dx * wavelength; % Шаг в метрах
            
            % Базовые амплитуды и фазы (без взаимодействия)
            base_matrix = zeros(12, 12, 2, 'like', 1i);
            
            % Вычисляем сдвиг фазы на столбец на основе theta_slider
            theta_shift_per_column = deg2rad(theta)*3.33; % Сдвиг фазы в радианах на один столбец
            
            % Число элементов в одной секции: 4 столбца * 12 рядов = 48
            elements_per_section = 48;
            
            % Секция 1: столбцы 1-4
            amplitude_variation = 1 + (amplitude_noise / 100) * rand(12, 4); % Случайный разброс амплитуды (0–50%)
            phase_base = (phase_noise / 100) * pi * rand(12, 4); % Случайная фаза ±(phase_noise/100)*π радиан
            for j = 1:4
                % Линейный сдвиг фазы для каждого столбца: (j-1) * theta_shift_per_column
                column_phase_shift = (j - 1) * theta_shift_per_column;
                if strcmp(app.ModeSwitch.Value, 'X-мода')
                    phase_1 = phase_base(:, j) + column_phase_shift; % Фаза первого диполя с учётом сдвига
                    phase_2 = phase_base(:, j) + column_phase_shift; % Фаза второго диполя (та же, что у первого)
                else % O-мода
                    phase_1 = phase_base(:, j) + column_phase_shift; % Фаза первого диполя с учётом сдвига
                    phase_2 = phase_base(:, j) + column_phase_shift + pi/2; % Фаза второго диполя (сдвиг на pi/2)
                end
                % Делим мощность на количество элементов в секции
                base_matrix(:, j, 1) = sqrt(powers(1) / elements_per_section / 4) * amplitude_factor * amplitude_variation(:, j) .* exp(1i * phase_1);
                base_matrix(:, j, 2) = sqrt(powers(1) / elements_per_section / 4) * amplitude_factor * amplitude_variation(:, j) .* exp(1i * phase_2);
            end
            
            % Секция 2: столбцы 5-8
            amplitude_variation = 1 + (amplitude_noise / 100) * rand(12, 4);
            phase_base = (phase_noise / 100) * pi * rand(12, 4);
            for j = 1:4
                % Линейный сдвиг фазы для столбцов 5-8: (j+3) * theta_shift_per_column
                column_phase_shift = (j + 3) * theta_shift_per_column;
                if strcmp(app.ModeSwitch.Value, 'X-мода')
                    phase_1 = phase_base(:, j) + column_phase_shift;
                    phase_2 = phase_base(:, j) + column_phase_shift;
                else % O-мода
                    phase_1 = phase_base(:, j) + column_phase_shift;
                    phase_2 = phase_base(:, j) + column_phase_shift + pi/2;
                end
                % Делим мощность на количество элементов в секции
                base_matrix(:, j+4, 1) = sqrt(powers(2) / elements_per_section / 4) * amplitude_factor * amplitude_variation(:, j) .* exp(1i * phase_1);
                base_matrix(:, j+4, 2) = sqrt(powers(2) / elements_per_section / 4) * amplitude_factor * amplitude_variation(:, j) .* exp(1i * phase_2);
            end
            
            % Секция 3: столбцы 9-12
            amplitude_variation = 1 + (amplitude_noise / 100) * rand(12, 4);
            phase_base = (phase_noise / 100) * pi * rand(12, 4);
            for j = 1:4
                column_phase_shift = (j + 7) * theta_shift_per_column;
                if strcmp(app.ModeSwitch.Value, 'X-мода')
                    phase_1 = phase_base(:, j) + column_phase_shift;
                    phase_2 = phase_base(:, j) + column_phase_shift;
                else % O-мода
                    phase_1 = phase_base(:, j) + column_phase_shift;
                    phase_2 = phase_base(:, j) + column_phase_shift + pi/2;
                end
                % Делим мощность на количество элементов в секции
                base_matrix(:, j+8, 1) = sqrt(powers(3) / elements_per_section / 4) * amplitude_factor * amplitude_variation(:, j) .* exp(1i * phase_1);
                base_matrix(:, j+8, 2) = sqrt(powers(3) / elements_per_section / 4) * amplitude_factor * amplitude_variation(:, j) .* exp(1i * phase_2);
            end
            
            % Присваиваем base_matrix в matrix (взаимодействие диполей удалено)
            matrix = base_matrix;
        end

        % Функция для вычисления амплитудной характеристики направленности
        function [psi_base, Fn_psi, theta_x, Fn_theta_x, theta_y, Fn_theta_y, theta_x_3d, theta_y_3d, Fn_3d, theta_sph, phi_sph, Fn_sph] = computeAmplitudePattern(app, frequency)
            % Параметры
            Nx_total = 12;
            Ny_total = 12;
            Nx = app.SelectedRow;
            Ny = app.SelectedColumn;
            d = 23.921;
            c = 3e8;
            ksi = 1;
            theta_slider = deg2rad(app.ThetaSlider.Value);
            theta_offset = 34;
        
            % Вычисление волнового числа
            f = frequency * 1e6;
            wavelength = c / f;
            k = 2 * pi / wavelength;
        
            % 1. Зависимость Fn от psi
            psi_0 = (Nx_total * k * d / 2) * cos(theta_slider + theta_offset - ksi);
            psi_base = linspace(-10*pi, 10*pi, 1000);
            psi_shifted = psi_base + psi_0;
            Fn_psi = abs(sin(psi_shifted) ./ (Nx_total * sin(psi_shifted / Nx_total)));
            Fn_psi(isnan(Fn_psi) | isinf(Fn_psi)) = 1;
        
            % 2. Зависимость Fn от theta_x и theta_y (для 3D-графика)
            theta_x_3d = linspace(-pi/2, pi/2, 361);
            theta_y_3d = linspace(-pi/2, pi/2, 361);
            [theta_x_3d, theta_y_3d] = meshgrid(theta_x_3d, theta_y_3d);
            Fn_3d = zeros(size(theta_x_3d));
        
            x_positions = d * ((1:Nx_total) - (Nx_total + 1)/2);
            y_positions = d * ((1:Ny_total) - (Ny_total + 1)/2);
        
            for i = 1:Nx_total
                for j = 1:Ny_total
                    for dipole_idx = 1:2
                        A_ij = abs(app.Matrix(i, j, dipole_idx));
                        phi_ij = angle(app.Matrix(i, j, dipole_idx));
                        dipole_angle = (dipole_idx == 1) * pi/4 + (dipole_idx == 2) * (pi/4 + pi/2);
                        phase = k * (x_positions(i) * cos(theta_x_3d + theta_offset - ksi) + y_positions(j) * cos(theta_y_3d + theta_offset - ksi));
                        directional_factor = cos(dipole_angle - atan2(cos(theta_y_3d), cos(theta_x_3d)));
                        Fn_3d = Fn_3d + A_ij * directional_factor .* exp(1i * (phi_ij + phase));
                    end
                end
            end
            Fn_3d = abs(Fn_3d);
            max_Fn_3d = max(Fn_3d(:));
            if max_Fn_3d > 0
                Fn_3d = Fn_3d / max_Fn_3d;
            else
                Fn_3d = zeros(size(Fn_3d));
            end
        
            % 2. Зависимость Fn от theta_x
            theta_x = linspace(-pi/2, pi/2, 361);
            nx_idx = round(1 + (Nx - 1) * (180 - 1) / (12 - 1));
            nx_idx = max(1, min(181, nx_idx));
            Fn_theta_x = Fn_3d(nx_idx, :);
            Fn_theta_x(isnan(Fn_theta_x)) = 0;
            max_Fn_theta_x = max(Fn_theta_x(:));
            if max_Fn_theta_x > 0
                Fn_theta_x = Fn_theta_x / max_Fn_theta_x;
            else
                Fn_theta_x = zeros(size(Fn_theta_x));
            end
            
            % 3. Зависимость Fn от theta_y
            theta_y = linspace(-pi/2, pi/2, 361);
            ny_idx = round(1 + (Ny - 1) * (180 - 1) / (12 - 1));
            ny_idx = max(1, min(181, ny_idx));
            Fn_theta_y = Fn_3d(:, ny_idx)';
            Fn_theta_y(isnan(Fn_theta_y)) = 0;
            max_Fn_theta_y = max(Fn_theta_y(:));
            if max_Fn_theta_y > 0
                Fn_theta_y = Fn_theta_y / max_Fn_theta_y;
            else
                Fn_theta_y = zeros(size(Fn_theta_y));
            end            
            % 4. Сферическая диаграмма направленности
            theta_x_3d = linspace(-pi/2, pi/2, 361);
            theta_y_3d = linspace(-pi/2, pi/2, 361);
            [theta_x_3d, theta_y_3d] = meshgrid(theta_x_3d, theta_y_3d);
            sin_theta_sq = sin(theta_x_3d).^2 + sin(theta_y_3d).^2;
            valid = sin_theta_sq <= 1;
            theta_sph = acos(sqrt(1 - sin_theta_sq .* valid));
            phi_sph = atan2(sin(theta_y_3d), sin(theta_x_3d));
            Fn_sph = Fn_3d;
            Fn_sph(~valid) = NaN;
            Fn_sph(isnan(Fn_sph)) = 0;
            max_Fn_sph = max(Fn_sph(:));
            if max_Fn_sph > 0
                Fn_sph = Fn_sph / max_Fn_sph;
            end
        end

        % Обновление графика "Амплитуда диполей по секциям"
        function updatePlot(app)
            % Вычисляем результирующую интенсивность с учётом фазы обоих диполей
            complex_sum = app.Matrix(:, :, 1) + app.Matrix(:, :, 2);
            intensity = abs(complex_sum).^2;
            intensity = flipud(intensity);
            
            if isempty(app.UIAxes) || ~isvalid(app.UIAxes)
                app.UIAxes = uiaxes(app.GridLayout);
                app.UIAxes.Layout.Row = 1;
                app.UIAxes.Layout.Column = 2;
            end
        
            cla(app.UIAxes);
            hh = image(app.UIAxes, intensity, 'CDataMapping', 'scaled');
            colormap(app.UIAxes, 'jet');
            colorbar(app.UIAxes);
            axis(app.UIAxes, [0.5 12.5 0.5 12.5]);
            
            xticks(app.UIAxes, 1:12);
            yticks(app.UIAxes, 1:12);
            yticklabels(app.UIAxes, string(12:-1:1));
            
            hold(app.UIAxes, 'on');
            for k = 0.5:1:13
                plot(app.UIAxes, [k k], [0 13], 'k-', 'LineWidth', 0.5);
                plot(app.UIAxes, [0 13], [k k], 'k-', 'LineWidth', 0.5);
            end
            hold(app.UIAxes, 'off');
            title(app.UIAxes, 'Температурная карта интенсивности сигналов с учётом фазы');
            
            set(hh, 'ButtonDownFcn', @(src, event) app.ElementClicked(src, event));
            
            drawnow;
        end

        % Callback для нажатия на элемент тепловой карты
        function ElementClicked(app, src, ~)
            clicked_point = app.UIAxes.CurrentPoint;
            x = round(clicked_point(1, 1));
            y = round(clicked_point(1, 2));
            
            if x >= 1 && x <= 12 && y >= 1 && y <= 12
                app.EditElementWindow(13 - y, x);
            end
        end

        % Создание окна редактирования элемента матрицы
        function EditElementWindow(app, i, j)
            editFig = uifigure('Name', ['Редактирование элемента (' num2str(i) ',' num2str(j) ')'], ...
                               'Position', [200 200 400 300], ...
                               'CloseRequestFcn', @(src, event) delete(src));
            
            editGrid = uigridlayout(editFig, [5 2]);
            editGrid.RowHeight = {30, 30, 30, 30, 30};
            editGrid.ColumnWidth = {'1x', '1x'};
            
            complex1 = app.Matrix(i, j, 1);
            complex2 = app.Matrix(i, j, 2);
            amp1 = abs(complex1);
            phase1 = angle(complex1) * 180 / pi;
            amp2 = abs(complex2);
            phase2 = angle(complex2) * 180 / pi;
            
            label1 = uilabel(editGrid, 'Text', 'Диполь 1: Амплитуда', 'HorizontalAlignment', 'right');
            label1.Layout.Row = 1;
            label1.Layout.Column = 1;
            
            label2 = uilabel(editGrid, 'Text', 'Диполь 1: Фаза (градусы)', 'HorizontalAlignment', 'right');
            label2.Layout.Row = 2;
            label2.Layout.Column = 1;
            
            ampField1 = uieditfield(editGrid, 'numeric', 'Value', amp1);
            ampField1.Layout.Row = 1;
            ampField1.Layout.Column = 2;
            
            phaseField1 = uieditfield(editGrid, 'numeric', 'Value', phase1);
            phaseField1.Layout.Row = 2;
            phaseField1.Layout.Column = 2;
            
            label3 = uilabel(editGrid, 'Text', 'Диполь 2: Амплитуда', 'HorizontalAlignment', 'right');
            label3.Layout.Row = 3;
            label3.Layout.Column = 1;
            
            label4 = uilabel(editGrid, 'Text', 'Диполь 2: Фаза (градусы)', 'HorizontalAlignment', 'right');
            label4.Layout.Row = 4;
            label4.Layout.Column = 1;
            
            ampField2 = uieditfield(editGrid, 'numeric', 'Value', amp2);
            ampField2.Layout.Row = 3;
            ampField2.Layout.Column = 2;
            
            phaseField2 = uieditfield(editGrid, 'numeric', 'Value', phase2);
            phaseField2.Layout.Row = 4;
            phaseField2.Layout.Column = 2;
            
            saveButton = uibutton(editGrid, 'Text', 'Сохранить');
            saveButton.Layout.Row = 5;
            saveButton.Layout.Column = [1 2];
            saveButton.ButtonPushedFcn = @(src, event) SaveElement(app, i, j, ampField1, phaseField1, ampField2, phaseField2, editFig);
        end

        % Сохранение изменений элемента матрицы
        function SaveElement(app, i, j, ampField1, phaseField1, ampField2, phaseField2, editFig)
            amp1 = ampField1.Value;
            phase1 = deg2rad(phaseField1.Value);
            amp2 = ampField2.Value;
            phase2 = deg2rad(phaseField2.Value);
            
            app.Matrix(i, j, 1) = amp1 * exp(1i * phase1);
            app.Matrix(i, j, 2) = amp2 * exp(1i * phase2);
            
            app.updatePlot();
            app.updateDirectionPlots();
            
            delete(editFig);
            
            app.appendMessage(['Элемент (' num2str(i) ',' num2str(j) ') обновлён']);
        end

        % Открытие окна редактирования всей матрицы через таблицу
        function OpenMatrixEditor(app)
            matrixFig = uifigure('Name', 'Редактирование матрицы', ...
                                 'Position', [200 200 800 600], ...
                                 'CloseRequestFcn', @(src, event) delete(src));
            
            matrixGrid = uigridlayout(matrixFig, [2 1]);
            matrixGrid.RowHeight = {'1x', 30};
            matrixGrid.ColumnWidth = {'1x'};
            
            numRows = 12 * 12;
            tableData = cell(numRows, 6);
            rowIdx = 1;
            for j = 1:12
                for i = 1:12
                    tableData{rowIdx, 1} = j;
                    tableData{rowIdx, 2} = i;
                    complex1 = app.Matrix(i, j, 1);
                    tableData{rowIdx, 3} = abs(complex1);
                    tableData{rowIdx, 4} = angle(complex1) * 180 / pi;
                    complex2 = app.Matrix(i, j, 2);
                    tableData{rowIdx, 5} = abs(complex2);
                    tableData{rowIdx, 6} = angle(complex2) * 180 / pi;
                    rowIdx = rowIdx + 1;
                end
            end
            
            matrixTable = uitable(matrixGrid);
            matrixTable.Data = tableData;
            matrixTable.ColumnName = {'Столбец', 'Ряд', 'Диполь 1: Амплитуда', 'Диполь 1: Фаза (град)', 'Диполь 2: Амплитуда', 'Диполь 2: Фаза (град)'};
            matrixTable.ColumnEditable = [false false true true true true];
            matrixTable.ColumnFormat = {'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric'};
            matrixTable.Layout.Row = 1;
            matrixTable.Layout.Column = 1;
            
            saveButton = uibutton(matrixGrid, 'Text', 'Сохранить');
            saveButton.Layout.Row = 2;
            saveButton.Layout.Column = 1;
            saveButton.ButtonPushedFcn = @(src, event) SaveMatrixFromTable(app, matrixTable, matrixFig);
        end

        % Сохранение изменений матрицы из таблицы
        function SaveMatrixFromTable(app, matrixTable, matrixFig)
            tableData = matrixTable.Data;
            
            for rowIdx = 1:size(tableData, 1)
                j = tableData{rowIdx, 1};
                i = tableData{rowIdx, 2};
                amp1 = tableData{rowIdx, 3};
                phase1 = deg2rad(tableData{rowIdx, 4});
                app.Matrix(i, j, 1) = amp1 * exp(1i * phase1);
                amp2 = tableData{rowIdx, 5};
                phase2 = deg2rad(tableData{rowIdx, 6});
                app.Matrix(i, j, 2) = amp2 * exp(1i * phase2);
            end
            
            app.updatePlot();
            app.updateDirectionPlots();
            
            delete(matrixFig);
            
            app.appendMessage('Матрица обновлена через таблицу');
        end

        % Обновление графиков направленности
        function updateDirectionPlots(app)
            frequency = app.FrequencySlider.Value;
            [psi_base, Fn_psi, theta_x, Fn_theta_x, theta_y, Fn_theta_y, theta_x_3d, theta_y_3d, Fn_3d, theta_sph, phi_sph, Fn_sph] = computeAmplitudePattern(app, frequency);
            
            cla(app.AxesPsi);
            plot(app.AxesPsi, psi_base, Fn_psi, 'b-');
            xlabel(app.AxesPsi, '\psi (радианы)');
            ylabel(app.AxesPsi, 'F_n(\psi)');
            title(app.AxesPsi, 'F_n( \psi)');
            grid(app.AxesPsi, 'on');
            
            cla(app.AxesThetaX);
            plot(app.AxesThetaX, rad2deg(theta_y), Fn_theta_y, 'r-');
            xlabel(app.AxesThetaX, '\theta_x (градусы)');
            ylabel(app.AxesThetaX, 'F_n(\theta_x) (нормализованная)');
            title(app.AxesThetaX, ['F_n(\theta_x) (ряд ' num2str(app.SelectedRow) ')']);
            grid(app.AxesThetaX, 'on');
            
            cla(app.AxesThetaY);
            plot(app.AxesThetaY, rad2deg(theta_x), Fn_theta_x, 'g-');
            xlabel(app.AxesThetaY, '\theta_y (градусы)');
            ylabel(app.AxesThetaY, 'F_n(\theta_y) (нормализованная)');
            title(app.AxesThetaY, ['F_n(\theta_y) (столбец ' num2str(app.SelectedColumn) ')']);
            grid(app.AxesThetaY, 'on');
            
            cla(app.AxesTheta3D);
            surf(app.AxesTheta3D, rad2deg(theta_y_3d), rad2deg(theta_x_3d), Fn_3d, 'EdgeColor', 'interp');
            xlabel(app.AxesTheta3D, '\theta_x (градусы)');
            ylabel(app.AxesTheta3D, '\theta_y (градусы)');
            zlabel(app.AxesTheta3D, 'F_n(\theta_x, \theta_y)');
            title(app.AxesTheta3D, 'F_n(\theta_x, \theta_y)');
            colormap(app.AxesTheta3D, 'jet');
            colorbar(app.AxesTheta3D);
            grid(app.AxesTheta3D, 'on');
            
            cla(app.AxesXY);
            r = Fn_sph;
            X = r .* sin(theta_sph) .* cos(phi_sph);
            Y = r .* sin(theta_sph) .* sin(phi_sph);
            Z = r .* cos(theta_sph);
            surf(app.AxesXY, X, Y, Z, r, 'EdgeColor', 'none');
            xlabel(app.AxesXY, 'Y');
            ylabel(app.AxesXY, 'X');
            zlabel(app.AxesXY, 'Z');
            title(app.AxesXY, 'Сферическая диаграмма направленности');
            colormap(app.AxesXY, 'jet');
            colorbar(app.AxesXY);
            axis(app.AxesXY, 'equal');
            grid(app.AxesXY, 'on');
            
            drawnow;
        end
        
        % Callback для слайдера мощности секции 1
        function PowerSlider1ValueChanged(app, ~)
            powers = [app.PowerSlider1.Value, app.PowerSlider2.Value, app.PowerSlider3.Value];
            frequency = app.FrequencySlider.Value;
            theta = app.ThetaSlider.Value;
            phase_noise = app.PhaseNoiseSlider.Value;
            amplitude_noise = app.AmplitudeNoiseSlider.Value;
            app.PowerLabel1.Text = ['Мощность ПКВ-250-1 = ' num2str(powers(1), '%.2f') ' кВт'];
            app.Matrix = initializeMatrix(app, powers, frequency, theta, phase_noise, amplitude_noise);
            updatePlot(app);
            app.appendMessage(['Изменена мощность ПКВ-250-1 на ' num2str(powers(1), '%.2f') ' кВт']);
            updateDirectionPlots(app);
        end

        % Callback для слайдера мощности секции 2
        function PowerSlider2ValueChanged(app, ~)
            powers = [app.PowerSlider1.Value, app.PowerSlider2.Value, app.PowerSlider3.Value];
            frequency = app.FrequencySlider.Value;
            theta = app.ThetaSlider.Value;
            phase_noise = app.PhaseNoiseSlider.Value;
            amplitude_noise = app.AmplitudeNoiseSlider.Value;
            app.PowerLabel2.Text = ['Мощность ПКВ-250-2 = ' num2str(powers(2), '%.2f') ' кВт'];
            app.Matrix = initializeMatrix(app, powers, frequency, theta, phase_noise, amplitude_noise);
            updatePlot(app);
            app.appendMessage(['Изменена мощность ПКВ-250-2 на ' num2str(powers(2), '%.2f') ' кВт']);
            updateDirectionPlots(app);
        end

        % Callback для слайдера мощности секции 3
        function PowerSlider3ValueChanged(app, ~)
            powers = [app.PowerSlider1.Value, app.PowerSlider2.Value, app.PowerSlider3.Value];
            frequency = app.FrequencySlider.Value;
            theta = app.ThetaSlider.Value;
            phase_noise = app.PhaseNoiseSlider.Value;
            amplitude_noise = app.AmplitudeNoiseSlider.Value;
            app.PowerLabel3.Text = ['Мощность ПКВ-250-3' num2str(powers(3), '%.2f') ' кВт'];
            app.Matrix = initializeMatrix(app, powers, frequency, theta, phase_noise, amplitude_noise);
            updatePlot(app);
            app.appendMessage(['Изменена мощность ПКВ-250-3']);
            updateDirectionPlots(app);
        end

        % Callback для слайдера частоты
        function FrequencySliderValueChanged(app, ~)
            powers = [app.PowerSlider1.Value, app.PowerSlider2.Value, app.PowerSlider3.Value];
            frequency = app.FrequencySlider.Value;
            theta = app.ThetaSlider.Value;
            phase_noise = app.PhaseNoiseSlider.Value;
            amplitude_noise = app.AmplitudeNoiseSlider.Value;
            app.FrequencyLabel.Text = ['Частота передатчиков = ' num2str(frequency, '%.2f') ' МГц'];
            app.Matrix = initializeMatrix(app, powers, frequency, theta, phase_noise, amplitude_noise);
            updatePlot(app);

            c = 3e8;
            f = frequency * 1e6;
            wavelength = c / f;
            k = 2 * pi / wavelength;

            app.appendMessage(['']);
            app.appendMessage(['Изменена частота на ' num2str(frequency, '%.2f') ' МГц']);
            app.appendMessage(['Длина волны ' num2str(wavelength) ' м']);
            app.appendMessage(['Волновое число k: ' num2str(k) ' радиан/м']);
            app.appendMessage(['']);

            updateDirectionPlots(app);
        end

        % Callback для слайдера угла theta
        function ThetaSliderValueChanged(app, ~)
            powers = [app.PowerSlider1.Value, app.PowerSlider2.Value, app.PowerSlider3.Value];
            frequency = app.FrequencySlider.Value;
            theta = app.ThetaSlider.Value;
            phase_noise = app.PhaseNoiseSlider.Value;
            amplitude_noise = app.AmplitudeNoiseSlider.Value;
            app.appendMessage(['']);
            app.ThetaLabel.Text = ['Угол сканирования = ' num2str(theta, '%.2f') ' градусов'];
            app.Matrix = initializeMatrix(app, powers, frequency, theta, phase_noise, amplitude_noise);
            updatePlot(app);
            app.appendMessage(['Изменён угол сканирования θ на значение ' num2str(theta, '%.2f') ' градусов']);
            app.appendMessage(['']);
            updateDirectionPlots(app);
        end

        % Callback для слайдера СКО фазы
        function PhaseNoiseSliderValueChanged(app, ~)
            powers = [app.PowerSlider1.Value, app.PowerSlider2.Value, app.PowerSlider3.Value];
            frequency = app.FrequencySlider.Value;
            theta = app.ThetaSlider.Value;
            phase_noise = app.PhaseNoiseSlider.Value;
            amplitude_noise = app.AmplitudeNoiseSlider.Value;
            app.PhaseNoiseLabel.Text = ['СКО фазы = ' num2str(phase_noise, '%.2f') ' %'];
            app.Matrix = initializeMatrix(app, powers, frequency, theta, phase_noise, amplitude_noise);
            updatePlot(app);
            app.appendMessage(['Изменено СКО фазы на ' num2str(phase_noise, '%.2f') ' %']);
            updateDirectionPlots(app);
        end

        % Callback для слайдера СКО амплитуды
        function AmplitudeNoiseSliderValueChanged(app, ~)
            powers = [app.PowerSlider1.Value, app.PowerSlider2.Value, app.PowerSlider3.Value];
            frequency = app.FrequencySlider.Value;
            theta = app.ThetaSlider.Value;
            phase_noise = app.PhaseNoiseSlider.Value;
            amplitude_noise = app.AmplitudeNoiseSlider.Value;
            app.AmplitudeNoiseLabel.Text = ['СКО амплитуды = ' num2str(amplitude_noise, '%.2f') ' %'];
            app.Matrix = initializeMatrix(app, powers, frequency, theta, phase_noise, amplitude_noise);
            updatePlot(app);
            app.appendMessage(['Изменено СКО амплитуды на ' num2str(amplitude_noise, '%.2f') ' %']);
            updateDirectionPlots(app);
        end

        % Callback для выпадающего списка ряда Nx
        function RowDropDownValueChanged(app, ~)
            app.SelectedRow = str2double(app.RowDropDown.Value);
            app.appendMessage(['']);
            app.appendMessage(['Выбран ряд Nx = ' num2str(app.SelectedRow)]);
            app.appendMessage(['']);
            updateDirectionPlots(app);
        end

        % Callback для выпадающего списка столбца Ny
        function ColumnDropDownValueChanged(app, ~)
            app.SelectedColumn = str2double(app.ColumnDropDown.Value);
            app.appendMessage(['']);
            app.appendMessage(['Выбран столбец Ny = ' num2str(app.SelectedColumn)]);
            updateDirectionPlots(app);
        end

        % Callback для кнопки отладки ЭМП
        function DirectionButtonPushed(app, ~)
            app.appendMessage('Запрос параметров ЭМП');
            frequency = app.FrequencySlider.Value;
            powers = [app.PowerSlider1.Value, app.PowerSlider2.Value, app.PowerSlider3.Value];
            
            c = 3e8;
            f = frequency * 1e6;
            wavelength = c / f;
            k = 2 * pi / wavelength;
            
            [~, ~, ~, ~, ~, ~, ~, ~, Fn_3d, ~, ~, ~] = computeAmplitudePattern(app, frequency);
            max_intensity_dB = 10 * log10(max(Fn_3d(:)));
            
            app.appendMessage(['Мощность ПКВ-250-1: ' num2str(powers(1)) ' кВт']);
            app.appendMessage(['Мощность ПКВ-250-2: ' num2str(powers(2)) ' кВт']);
            app.appendMessage(['Мощность ПКВ-250-3: ' num2str(powers(3)) ' кВт']);
            app.appendMessage(['']);
            app.appendMessage(['Частота излучения: ' num2str(frequency) ' МГц']);
            app.appendMessage(['Длина волны: ' num2str(wavelength) ' м']);
            app.appendMessage(['Волновое число k: ' num2str(k) ' радиан/м']);
            app.appendMessage(['']);
            app.appendMessage(['Максимум интенсивности в дальней зоне: ' num2str(max_intensity_dB) ' дБ']);
            app.appendMessage(['']);
            drawnow;
        end
        
        % Callback для тумблера режима (X-мода или O-мода)
        function ModeSwitchValueChanged(app, ~)
            powers = [app.PowerSlider1.Value, app.PowerSlider2.Value, app.PowerSlider3.Value];
            frequency = app.FrequencySlider.Value;
            theta = app.ThetaSlider.Value;
            phase_noise = app.PhaseNoiseSlider.Value;
            amplitude_noise = app.AmplitudeNoiseSlider.Value;
            app.appendMessage(['Изменён режим на ' app.ModeSwitch.Value]);
            app.Matrix = initializeMatrix(app, powers, frequency, theta, phase_noise, amplitude_noise);
            updatePlot(app);
            updateDirectionPlots(app);
        end
        
        % Callback для кнопки построения диаграммы направленности
        function DNPButtonPushed(app, ~)
            app.appendMessage('Построение характеристик направленности...');
            updateDirectionPlots(app);
        end
    end

    methods (Access = public)

        % Конструктор приложения
        function app = DipoleMatrixApp14
            app.UIFigure = uifigure('Name', 'Модель нагревного стенда СУРА');
            app.UIFigure.Position = [0 0 1680 920];

            app.GridLayout = uigridlayout(app.UIFigure, [2 3]);
            app.GridLayout.RowHeight = {'1x', '1x'};
            app.GridLayout.ColumnWidth = {'1x', '2x', '5x'};
            app.GridLayout.RowSpacing = 5;
            app.GridLayout.ColumnSpacing = 5;
            app.GridLayout.BackgroundColor = [0.7 0.7 0.7];

            controlsPanel = uipanel(app.GridLayout);
            controlsPanel.Layout.Row = 1;
            controlsPanel.Layout.Column = 1;
            controlsPanel.Title = 'Блок управления стендом СУРА';
            controlsPanel.BackgroundColor = [0.95 0.95 0.95];

            app.PowerLabel1 = uilabel(controlsPanel);
            app.PowerLabel1.Position = [20 420 200 22];
            app.PowerLabel1.Text = 'Мощность ПКВ-250-1, кВт';
            
            app.PowerSlider1 = uislider(controlsPanel);
            app.PowerSlider1.Limits = [0 300];
            app.PowerSlider1.Value = 6.25;
            app.PowerSlider1.Position = [20 415 200 3];
            app.PowerSlider1.MajorTicks = [0:50:300];
            app.PowerSlider1.ValueChangedFcn = @(src, event) app.PowerSlider1ValueChanged;

            app.PowerLabel2 = uilabel(controlsPanel);
            app.PowerLabel2.Position = [20 365 200 22];
            app.PowerLabel2.Text = 'Мощность ПКВ-250-2, кВт';
            
            app.PowerSlider2 = uislider(controlsPanel);
            app.PowerSlider2.Limits = [0 300];
            app.PowerSlider2.Value = 6.25;
            app.PowerSlider2.Position = [20 360 200 3];
            app.PowerSlider2.MajorTicks = [0:50:300];
            app.PowerSlider2.ValueChangedFcn = @(src, event) app.PowerSlider2ValueChanged;

            app.PowerLabel3 = uilabel(controlsPanel);
            app.PowerLabel3.Position = [20 310 200 22];
            app.PowerLabel3.Text = 'Мощность ПКВ-250-3, кВт';
            
            app.PowerSlider3 = uislider(controlsPanel);
            app.PowerSlider3.Limits = [0 300];
            app.PowerSlider3.Value = 6.25;
            app.PowerSlider3.Position = [20 305 200 3];
            app.PowerSlider3.MajorTicks = [0:50:300];
            app.PowerSlider3.ValueChangedFcn = @(src, event) app.PowerSlider3ValueChanged;

            app.FrequencyLabel = uilabel(controlsPanel);
            app.FrequencyLabel.Position = [20 255 200 22];
            app.FrequencyLabel.Text = 'Частота передатчиков, МГц';
            
            app.FrequencySlider = uislider(controlsPanel);
            app.FrequencySlider.Limits = [3 10];
            app.FrequencySlider.Value = 5;
            app.FrequencySlider.Position = [20 250 200 3];
            app.FrequencySlider.MajorTicks = [3:1:10];
            app.FrequencySlider.ValueChangedFcn = @(src, event) app.FrequencySliderValueChanged;
            
            app.ThetaLabel = uilabel(controlsPanel);
            app.ThetaLabel.Position = [20 200 200 22];
            app.ThetaLabel.Text = 'Угол луча, градусы';
            
            app.ThetaSlider = uislider(controlsPanel);
            app.ThetaSlider.Limits = [-40 40];
            app.ThetaSlider.Value = 0;
            app.ThetaSlider.Position = [20 195 200 3];
            app.ThetaSlider.MajorTicks = [-40:10:40];
            app.ThetaSlider.ValueChangedFcn = @(src, event) app.ThetaSliderValueChanged;
            
            app.RowLabel = uilabel(controlsPanel);
            app.RowLabel.Position = [8 140 200 22];
            app.RowLabel.Text = 'Выбор ряда Nx';
            
            app.RowDropDown = uidropdown(controlsPanel);
            app.RowDropDown.Items = string(1:12);
            app.RowDropDown.Value = '1';
            app.RowDropDown.Position = [8 120 90 22];
            app.RowDropDown.ValueChangedFcn = @(src, event) app.RowDropDownValueChanged;
            
            app.ColumnLabel = uilabel(controlsPanel);
            app.ColumnLabel.Position = [128 140 200 22];
            app.ColumnLabel.Text = 'Выбор столбца Ny';
            
            app.ColumnDropDown = uidropdown(controlsPanel);
            app.ColumnDropDown.Items = string(1:12);
            app.ColumnDropDown.Value = '1';
            app.ColumnDropDown.Position = [128 120 90 22];
            app.ColumnDropDown.ValueChangedFcn = @(src, event) app.ColumnDropDownValueChanged;
            
            app.ModeSwitch = uiswitch(controlsPanel, 'slider');
            app.ModeSwitch.Position = [88 95 90 22];
            app.ModeSwitch.Items = {'X-мода', 'O-мода'};
            app.ModeSwitch.Value = 'X-мода';
            app.ModeSwitch.ValueChangedFcn = @(src, event) app.ModeSwitchValueChanged;

            % Слайдер и метка для СКО фазы
            app.PhaseNoiseLabel = uilabel(controlsPanel);
            app.PhaseNoiseLabel.Position = [20 75 200 22];
            app.PhaseNoiseLabel.Text = 'СКО фазы, %';
            
            app.PhaseNoiseSlider = uislider(controlsPanel);
            app.PhaseNoiseSlider.Limits = [0 50];
            app.PhaseNoiseSlider.Value = 0;
            app.PhaseNoiseSlider.Position = [20 70 200 3];
            app.PhaseNoiseSlider.MajorTicks = [0:10:50];
            app.PhaseNoiseSlider.ValueChangedFcn = @(src, event) app.PhaseNoiseSliderValueChanged;

            % Слайдер и метка для СКО амплитуды
            app.AmplitudeNoiseLabel = uilabel(controlsPanel);
            app.AmplitudeNoiseLabel.Position = [20 20 200 22];
            app.AmplitudeNoiseLabel.Text = 'СКО амплитуды, %';
            
            app.AmplitudeNoiseSlider = uislider(controlsPanel);
            app.AmplitudeNoiseSlider.Limits = [0 50];
            app.AmplitudeNoiseSlider.Value = 0;
            app.AmplitudeNoiseSlider.Position = [20 15 200 3];
            app.AmplitudeNoiseSlider.MajorTicks = [0:10:50];
            app.AmplitudeNoiseSlider.ValueChangedFcn = @(src, event) app.AmplitudeNoiseSliderValueChanged;

            app.DirectionButton = uicontrol(controlsPanel, 'Style', 'pushbutton');
            app.DirectionButton.String = 'Отладочный вывод ЭМП';
            app.DirectionButton.Position = [20 -30 90 22];
            app.DirectionButton.Callback = @(src, event) app.DirectionButtonPushed(app);
            


            intensityPanel = uipanel(app.GridLayout);
            intensityPanel.Layout.Row = 1;
            intensityPanel.Layout.Column = 2;
            intensityPanel.Title = 'Блок управления и мониторинга параметров сигналов излучателей';
            intensityPanel.BackgroundColor = [0.95 0.95 0.95];

            app.UIAxes = uiaxes(intensityPanel);
            app.UIAxes.Position = [10 40 400 270];
            
            dataButton = uibutton(intensityPanel, 'Text', 'Ввод данных');
            dataButton.Position = [25 140 100 22];
            dataButton.ButtonPushedFcn = @(src, event) app.OpenMatrixEditor();

            debugPanel = uipanel(app.GridLayout);
            debugPanel.Layout.Row = 2;
            debugPanel.Layout.Column = [1 2];
            debugPanel.Title = 'Центр контроля и отладки';
            debugPanel.BackgroundColor = [0.95 0.95 0.95];
            
            app.TextArea = uitextarea(debugPanel);
            app.TextArea.Position = [10 10 580 150];
            app.TextArea.Editable = 'off';

            directionPanel = uipanel(app.GridLayout);
            directionPanel.Layout.Row = [1 2];
            directionPanel.Layout.Column = 3;
            directionPanel.Title = 'Поле в дальней зоне';
            directionPanel.BackgroundColor = [0.95 0.95 0.95];
            
            directionGrid = uigridlayout(directionPanel, [2 3]);
            directionGrid.RowHeight = {'1x', '1x'};
            directionGrid.ColumnWidth = {'1x', '1x', '1x'};
            directionGrid.RowSpacing = 5;
            directionGrid.ColumnSpacing = 5;
            
            app.AxesPsi = uiaxes(directionGrid);
            app.AxesPsi.Layout.Row = 1;
            app.AxesPsi.Layout.Column = 1;
            
            app.AxesThetaX = uiaxes(directionGrid);
            app.AxesThetaX.Layout.Row = 1;
            app.AxesThetaX.Layout.Column = 2;
            
            app.AxesThetaY = uiaxes(directionGrid);
            app.AxesThetaY.Layout.Row = 1;
            app.AxesThetaY.Layout.Column = 3;
            
            app.AxesTheta3D = uiaxes(directionGrid);
            app.AxesTheta3D.Layout.Row = 2;
            app.AxesTheta3D.Layout.Column = [1 2];
            
            app.AxesXY = uiaxes(directionGrid);
            app.AxesXY.Layout.Row = 2;
            app.AxesXY.Layout.Column = 3;

            startupFcn(app);
        end
    end

    methods (Static)
        function run()
            app = DipoleMatrixApp14;
            disp('Интерфейс запущен.');
        end
    end
end