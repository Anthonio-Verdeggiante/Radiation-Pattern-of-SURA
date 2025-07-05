classdef DipoleMatrixApp12 < matlab.apps.AppBase

    % Свойства приложения
    properties (Access = public)
        UIFigure          matlab.ui.Figure
        GridLayout        matlab.ui.container.GridLayout % Сетка для разделения интерфейса
        PowerSlider1      matlab.ui.control.Slider % Мощность секции 1
        PowerSlider2      matlab.ui.control.Slider % Мощность секции 2
        PowerSlider3      matlab.ui.control.Slider % Мощность секции 3
        FrequencySlider   matlab.ui.control.Slider % Единый слайдер частоты
        ThetaSlider       matlab.ui.control.Slider % Слайдер для угла theta
        RowDropDown       matlab.ui.control.DropDown % Выпадающий список для выбора ряда Nx
        ColumnDropDown    matlab.ui.control.DropDown % Выпадающий список для выбора столбца Ny
        PowerLabel1       matlab.ui.control.Label
        PowerLabel2       matlab.ui.control.Label
        PowerLabel3       matlab.ui.control.Label
        FrequencyLabel    matlab.ui.control.Label
        ThetaLabel        matlab.ui.control.Label  % Метка для угла theta
        RowLabel          matlab.ui.control.Label  % Метка для выбора ряда Nx
        ColumnLabel       matlab.ui.control.Label  % Метка для выбора столбца Ny
        ModeSwitch        matlab.ui.control.Switch % Тумблер для выбора режима (X-мода или O-мода)
        UIAxes            matlab.ui.control.UIAxes % Для интенсивности матрицы
        ScaleAxes         matlab.ui.control.UIAxes % Ось для цветовой шкалы (новое свойство)
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
            app.Matrix = initializeMatrix(app, [300, 300, 300], app.FrequencySlider.Value, app.ThetaSlider.Value);
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
        function matrix = initializeMatrix(app, powers, frequency, theta)
            matrix = zeros(12, 12, 2); % Матрица 12x12x2 для двух диполей на элемент
            c = 3e8; % Скорость света, м/с
            dipole_length = 33.830; % Длина диполя, м (изменено с 25 м на 34 м)
            resonant_frequency = c / (2 * dipole_length); % Резонансная частота
            resonant_frequency_MHz = resonant_frequency / 1e6; % Перевод в МГц (~4.41 МГц)
            
            % Резонансная характеристика (гауссова функция)
            bandwidth = 5.0; % Полоса пропускания, МГц
            amplitude_factor = exp(-((frequency - resonant_frequency_MHz) / bandwidth).^2);
            
            % Шаг между диполями в метрах
            d = 23.921; % Шаг в долях длины волны
            wavelength = c / (frequency * 1e6);
            
            % Базовые амплитуды и фазы
            base_matrix = zeros(12, 12, 2, 'like', 1i);
            
            % Вычисляем сдвиг фазы на столбец на основе theta_slider
            theta_shift_per_row = deg2rad(theta)*3.33; % Сдвиг фазы в радианах на один столбец
            
            % Число перекрестных диполей в каждой секции (12 строк × 4 столбца = 48)
            dipoles_per_section = 48;
            elements_per_dipole = 2; % Перекрестный диполь состоит из 2 элементарных диполей
            
            % Секция 1: столбцы 1-4
            amplitude_variation = 1 + 0 * rand(12, 4);
            phase_base = 2 * pi * frequency * (dipole_length / c) * 1*ones(12, 4);
            for i = 1:4
                row_phase_shift = (i - 1) * theta_shift_per_row;
                if strcmp(app.ModeSwitch.Value, 'X-мода')
                    phase_1 = phase_base(:, i) + row_phase_shift;
                    phase_2 = phase_base(:, i) + row_phase_shift;
                else % O-мода
                    phase_1 = phase_base(:, i) + row_phase_shift;
                    phase_2 = phase_base(:, i) + row_phase_shift + pi/2;
                end
                % Корректная амплитуда на элементарный диполь
                dipole_power = powers(1) / (dipoles_per_section * elements_per_dipole); % Мощность на элементарный диполь
                amplitude = sqrt(dipole_power); % Амплитуда пропорциональна корню из мощности
                base_matrix(:, i, 1) = amplitude * amplitude_factor * amplitude_variation(:, i) .* exp(1i * phase_1);
                base_matrix(:, i, 2) = amplitude * amplitude_factor * amplitude_variation(:, i) .* exp(1i * phase_2);
            end
            
            % Секция 2: столбцы 5-8
            amplitude_variation = 1 + 0 * rand(12, 4);
            phase_base = 2 * pi * frequency * (dipole_length / c) * 1*ones(12, 4);
            for i = 1:4
                row_phase_shift = (i + 3) * theta_shift_per_row;
                if strcmp(app.ModeSwitch.Value, 'X-мода')
                    phase_1 = phase_base(:, i) + row_phase_shift;
                    phase_2 = phase_base(:, i) + row_phase_shift;
                else % O-мода
                    phase_1 = phase_base(:, i) + row_phase_shift;
                    phase_2 = phase_base(:, i) + row_phase_shift + pi/2;
                end
                dipole_power = powers(2) / (dipoles_per_section * elements_per_dipole); % Мощность на элементарный диполь
                amplitude = sqrt(dipole_power); % Амплитуда пропорциональна корню из мощности
                base_matrix(:, i+4, 1) = amplitude * amplitude_factor * amplitude_variation(:, i) .* exp(1i * phase_1);
                base_matrix(:, i+4, 2) = amplitude * amplitude_factor * amplitude_variation(:, i) .* exp(1i * phase_2);
            end
            
            % Секция 3: столбцы 9-12
            amplitude_variation = 1 + 0 * rand(12, 4);
            phase_base = 2 * pi * frequency * (dipole_length / c) * 1*ones(12, 4);
            for i = 1:4
                row_phase_shift = (i + 7) * theta_shift_per_row;
                if strcmp(app.ModeSwitch.Value, 'X-мода')
                    phase_1 = phase_base(:, i) + row_phase_shift;
                    phase_2 = phase_base(:, i) + row_phase_shift;
                else % O-мода
                    phase_1 = phase_base(:, i) + row_phase_shift;
                    phase_2 = phase_base(:, i) + row_phase_shift + pi/2;
                end
                dipole_power = powers(3) / (dipoles_per_section * elements_per_dipole); % Мощность на элементарный диполь
                amplitude = sqrt(dipole_power); % Амплитуда пропорциональна корню из мощности
                base_matrix(:, i+8, 1) = amplitude * amplitude_factor * amplitude_variation(:, i) .* exp(1i * phase_1);
                base_matrix(:, i+8, 2) = amplitude * amplitude_factor * amplitude_variation(:, i) .* exp(1i * phase_2);
            end
            
            % Устанавливаем итоговую матрицу без учета взаимодействия
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
            nx_idx = round(1 + (Nx - 1) * (180 - 1) / (12 - 1)); % 181 - 1 = 180
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
            intensity = flipud(intensity); % Отражаем матрицу по вертикали для корректного отображения
            
            if isempty(app.UIAxes) || ~isvalid(app.UIAxes)
                app.UIAxes = uiaxes(app.GridLayout);
                app.UIAxes.Layout.Row = 1;
                app.UIAxes.Layout.Column = 2;
            end
        
            % Очищаем оси
            cla(app.UIAxes);
            % Создаём тепловую карту с помощью image
            hh = image(app.UIAxes, intensity, 'CDataMapping', 'scaled');
            colormap(app.UIAxes, 'jet');
        
            % Устанавливаем абсолютные границы цветовой шкалы
            app.UIAxes.CLim = [0 12.18];
            axis(app.UIAxes, [0.5 12.5 0.5 12.5]); % Границы для 12x12 матрицы
            
            % Настраиваем метки осей
            xticks(app.UIAxes, 1:12); % Метки по горизонтальной оси: 1, 2, ..., 12
            yticks(app.UIAxes, 1:12); % Метки по вертикальной оси: 1, 2, ..., 12 снизу вверх
            yticklabels(app.UIAxes, string(12:-1:1)); % Меняем метки по вертикали: 12 сверху, 1 снизу
            
            % Добавляем сетку
            hold(app.UIAxes, 'on');
            for k = 0.5:1:13 % Обновляем диапазон сетки для 12 элементов
                plot(app.UIAxes, [k k], [0 13], 'k-', 'LineWidth', 0.5);
                plot(app.UIAxes, [0 13], [k k], 'k-', 'LineWidth', 0.5);
            end
            hold(app.UIAxes, 'off');
            title(app.UIAxes, 'Температурная карта интенсивности сигналов с учётом фазы');
            
            % Добавляем обработчик события нажатия на тепловую карту
            set(hh, 'ButtonDownFcn', @(src, event) app.ElementClicked(src, event));
            
            % Добавляем цветовую шкалу как отдельный градиент
            % Создаём оси для шкалы справа от графика
            scaleAxes = axes('Parent', app.UIFigure, 'Position', [0.75 0.1 0.02 0.4]); % Позиция: [x y ширина высота]
            colormap(scaleAxes, 'jet');
            % Создаём градиент от 0 до 1
            gradient = (0:0.01:1)';
            image(scaleAxes, 1, gradient, (1:size(gradient, 1))', 'CDataMapping', 'scaled');
            set(scaleAxes, 'YDir', 'normal', 'XTick', [], 'YTick', [0 0.25 0.5 0.75 1], ...
                           'YTickLabel', {'0', '~0.885', '~1.77', '~2.655', '3.54'}, ...
                           'TickLabelInterpreter', 'latex');
            ylabel(scaleAxes, 'Мощность (кВт)', 'Rotation', 0, 'HorizontalAlignment', 'right');
            axis(scaleAxes, 'off'); % Убираем оси
            
            drawnow;
        end

        % Callback для нажатия на элемент тепловой карты
        function ElementClicked(app, src, ~)
            % Получаем координаты нажатия
            clicked_point = app.UIAxes.CurrentPoint;
            x = round(clicked_point(1, 1));
            y = round(clicked_point(1, 2));
            
            % Проверяем, что нажатие произошло внутри допустимых границ матрицы
            if x >= 1 && x <= 12 && y >= 1 && y <= 12
                % Открываем окно редактирования для элемента (y, x)
                app.EditElementWindow(13 - y, x);
            end
        end

        % Создание окна редактирования элемента матрицы
        function EditElementWindow(app, i, j)
            % Создаём новое окно
            editFig = uifigure('Name', ['Редактирование элемента (' num2str(i) ',' num2str(j) ')'], ...
                               'Position', [200 200 400 300], ...
                               'CloseRequestFcn', @(src, event) delete(src));
            
            % Создаём сетку для размещения элементов интерфейса
            editGrid = uigridlayout(editFig, [5 2]);
            editGrid.RowHeight = {30, 30, 30, 30, 30};
            editGrid.ColumnWidth = {'1x', '1x'};
            
            % Текущие значения комплексных чисел
            complex1 = app.Matrix(i, j, 1);
            complex2 = app.Matrix(i, j, 2);
            amp1 = abs(complex1);
            phase1 = angle(complex1) * 180 / pi; % Переводим фазу в градусы
            amp2 = abs(complex2);
            phase2 = angle(complex2) * 180 / pi; % Переводим фазу в градусы
            
            % Поля для первого диполя
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
            
            % Поля для второго диполя
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
            
            % Кнопка сохранения
            saveButton = uibutton(editGrid, 'Text', 'Сохранить');
            saveButton.Layout.Row = 5;
            saveButton.Layout.Column = [1 2];
            saveButton.ButtonPushedFcn = @(src, event) SaveElement(app, i, j, ampField1, phaseField1, ampField2, phaseField2, editFig);
        end

        % Сохранение изменений элемента матрицы
        function SaveElement(app, i, j, ampField1, phaseField1, ampField2, phaseField2, editFig)
            % Получаем новые значения
            amp1 = ampField1.Value;
            phase1 = deg2rad(phaseField1.Value); % Переводим фазу из градусов в радианы
            amp2 = ampField2.Value;
            phase2 = deg2rad(phaseField2.Value); % Переводим фазу из градусов в радианы
            
            % Обновляем комплексные числа
            app.Matrix(i, j, 1) = amp1 * exp(1i * phase1);
            app.Matrix(i, j, 2) = amp2 * exp(1i * phase2);
            
            % Обновляем графики
            app.updatePlot();
            app.updateDirectionPlots();
            
            % Закрываем окно редактирования
            delete(editFig);
            
            % Сообщаем об успешном изменении
            app.appendMessage(['Элемент (' num2str(i) ',' num2str(j) ') обновлён']);
        end

        % Открытие окна редактирования всей матрицы через таблицу
        function OpenMatrixEditor(app)
            % Создаём новое окно
            matrixFig = uifigure('Name', 'Редактирование матрицы', ...
                                 'Position', [200 200 800 600], ...
                                 'CloseRequestFcn', @(src, event) delete(src));
            
            % Создаём сетку для размещения таблицы и кнопки
            matrixGrid = uigridlayout(matrixFig, [2 1]);
            matrixGrid.RowHeight = {'1x', 30};
            matrixGrid.ColumnWidth = {'1x'};
            
            % Подготовка данных для таблицы
            numRows = 12 * 12; % 144 элемента
            tableData = cell(numRows, 6);
            rowIdx = 1;
            for j = 1:12 % Сначала перебираем столбцы (слева направо)
                for i = 1:12 % Затем перебираем строки (сверху вниз)
                    % Столбец и ряд
                    tableData{rowIdx, 1} = j; % Столбец
                    tableData{rowIdx, 2} = i; % Ряд
                    % Диполь 1
                    complex1 = app.Matrix(i, j, 1);
                    tableData{rowIdx, 3} = abs(complex1); % Амплитуда
                    tableData{rowIdx, 4} = angle(complex1) * 180 / pi; % Фаза в градусах
                    % Диполь 2
                    complex2 = app.Matrix(i, j, 2);
                    tableData{rowIdx, 5} = abs(complex2); % Амплитуда
                    tableData{rowIdx, 6} = angle(complex2) * 180 / pi; % Фаза в градусах
                    rowIdx = rowIdx + 1;
                end
            end
            
            % Создаём таблицу
            matrixTable = uitable(matrixGrid);
            matrixTable.Data = tableData;
            matrixTable.ColumnName = {'Столбец', 'Ряд', 'Диполь 1: Амплитуда', 'Диполь 1: Фаза (град)', 'Диполь 2: Амплитуда', 'Диполь 2: Фаза (град)'};
            matrixTable.ColumnEditable = [false false true true true true]; % Столбец и Ряд нельзя редактировать
            matrixTable.ColumnFormat = {'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric'};
            matrixTable.Layout.Row = 1;
            matrixTable.Layout.Column = 1;
            
            % Кнопка сохранения
            saveButton = uibutton(matrixGrid, 'Text', 'Сохранить');
            saveButton.Layout.Row = 2;
            saveButton.Layout.Column = 1;
            saveButton.ButtonPushedFcn = @(src, event) SaveMatrixFromTable(app, matrixTable, matrixFig);
        end

        % Сохранение изменений матрицы из таблицы
        function SaveMatrixFromTable(app, matrixTable, matrixFig)
            % Получаем данные из таблицы
            tableData = matrixTable.Data;
            
            % Обновляем матрицу app.Matrix
            for rowIdx = 1:size(tableData, 1)
                j = tableData{rowIdx, 1}; % Столбец
                i = tableData{rowIdx, 2}; % Ряд
                % Диполь 1
                amp1 = tableData{rowIdx, 3};
                phase1 = deg2rad(tableData{rowIdx, 4}); % Фаза в радианах
                app.Matrix(i, j, 1) = amp1 * exp(1i * phase1);
                % Диполь 2
                amp2 = tableData{rowIdx, 5};
                phase2 = deg2rad(tableData{rowIdx, 6}); % Фаза в радианах
                app.Matrix(i, j, 2) = amp2 * exp(1i * phase2);
            end
            
            % Обновляем графики
            app.updatePlot();
            app.updateDirectionPlots();
            
            % Закрываем окно редактирования
            delete(matrixFig);
            
            % Сообщаем об успешном изменении
            app.appendMessage('Матрица обновлена через таблицу');
        end

        % Обновление графиков направленности
        function updateDirectionPlots(app)
            frequency = app.FrequencySlider.Value;
            [psi_base, Fn_psi, theta_x, Fn_theta_x, theta_y, Fn_theta_y, theta_x_3d, theta_y_3d, Fn_3d, theta_sph, phi_sph, Fn_sph] = computeAmplitudePattern(app, frequency);
            
            % График 1: Fn(psi)
            cla(app.AxesPsi);
            plot(app.AxesPsi, psi_base, Fn_psi, 'b-');
            xlabel(app.AxesPsi, '\psi (радианы)');
            ylabel(app.AxesPsi, 'F_n(\psi)');
            title(app.AxesPsi, 'F_n( \psi)');
            grid(app.AxesPsi, 'on');
            
            % График 2: Fn(theta_x)
            cla(app.AxesThetaX);
            plot(app.AxesThetaX, rad2deg(theta_y), Fn_theta_y, 'r-');
            xlabel(app.AxesThetaX, '\theta_x (градусы)');
            ylabel(app.AxesThetaX, 'F_n(\theta_x) (нормализованная)');
            title(app.AxesThetaX, ['F_n(\theta_x) (ряд ' num2str(app.SelectedRow) ')']);
            grid(app.AxesThetaX, 'on');
            
            % График 3: Fn(theta_y)
            cla(app.AxesThetaY);
            plot(app.AxesThetaY, rad2deg(theta_x), Fn_theta_x, 'g-');
            xlabel(app.AxesThetaY, '\theta_y (градусы)');
            ylabel(app.AxesThetaY, 'F_n(\theta_y) (нормализованная)');
            title(app.AxesThetaY, ['F_n(\theta_y) (столбец ' num2str(app.SelectedColumn) ')']);
            grid(app.AxesThetaY, 'on');
            
            % График 4: Fn(theta_x, theta_y)
            cla(app.AxesTheta3D);
            surf(app.AxesTheta3D, rad2deg(theta_y_3d), rad2deg(theta_x_3d), Fn_3d, 'EdgeColor', 'interp');
            xlabel(app.AxesTheta3D, '\theta_x (градусы)');
            ylabel(app.AxesTheta3D, '\theta_y (градусы)');
            zlabel(app.AxesTheta3D, 'F_n(\theta_x, \theta_y)');
            title(app.AxesTheta3D, 'F_n(\theta_x, \theta_y)');
            colormap(app.AxesTheta3D, 'jet');
            colorbar(app.AxesTheta3D);
            grid(app.AxesTheta3D, 'on');
            
            % График 5: Сферическая диаграмма направленности
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
            app.PowerLabel1.Text = ['Мощность ПКВ-250-1 = ' num2str(powers(1), '%.2f') ' кВт'];
            app.Matrix = initializeMatrix(app, powers, frequency, theta);
            updatePlot(app);
            app.appendMessage(['Изменена мощность ПКВ-250-1 на ' num2str(powers(1), '%.2f') ' кВт']);
            updateDirectionPlots(app); % Обновляем графики направленности
        end

        % Callback для слайдера мощности секции 2
        function PowerSlider2ValueChanged(app, ~)
            powers = [app.PowerSlider1.Value, app.PowerSlider2.Value, app.PowerSlider3.Value];
            frequency = app.FrequencySlider.Value;
            theta = app.ThetaSlider.Value;
            app.PowerLabel2.Text = ['Мощность ПКВ-250-2 = ' num2str(powers(2), '%.2f') ' кВт'];
            app.Matrix = initializeMatrix(app, powers, frequency, theta);
            updatePlot(app);
            app.appendMessage(['Изменена мощность ПКВ-250-2 на ' num2str(powers(2), '%.2f') ' кВт']);
            updateDirectionPlots(app); % Обновляем графики направленности
        end

        % Callback для слайдера мощности секции 3
        function PowerSlider3ValueChanged(app, ~)
            powers = [app.PowerSlider1.Value, app.PowerSlider2.Value, app.PowerSlider3.Value];
            frequency = app.FrequencySlider.Value;
            theta = app.ThetaSlider.Value;
            app.PowerLabel3.Text = ['Мощность ПКВ-250-3 = ' num2str(powers(3), '%.2f') ' кВт'];
            app.Matrix = initializeMatrix(app, powers, frequency, theta);
            updatePlot(app);
            app.appendMessage(['Изменена мощность ПКВ-250-3 на ' num2str(powers(3), '%.2f') ' кВт']);
            updateDirectionPlots(app); % Обновляем графики направленности
        end

        % Callback для слайдера частоты
        function FrequencySliderValueChanged(app, ~)
            powers = [app.PowerSlider1.Value, app.PowerSlider2.Value, app.PowerSlider3.Value];
            frequency = app.FrequencySlider.Value;
            theta = app.ThetaSlider.Value;
            app.FrequencyLabel.Text = ['Частота передатчиков = ' num2str(frequency, '%.2f') ' МГц'];
            app.Matrix = initializeMatrix(app, powers, frequency, theta);
            updatePlot(app);

            % Вычисляем длину волны и волновое число
            c = 3e8; % Скорость света, м/с
            f = frequency * 1e6; % Частота в Гц (перевод из МГц)
            wavelength = c / f; % Длина волны, м
            k = 2 * pi / wavelength; % Волновое число, радиан/м

            % Выводим сообщения
            app.appendMessage(['']);
            app.appendMessage(['Изменена частота на ' num2str(frequency, '%.2f') ' МГц']);
            app.appendMessage(['Длина волны: ' num2str(wavelength) ' м']);
            app.appendMessage(['Волновое число k: ' num2str(k) ' радиан/м']);
            app.appendMessage(['']);

            updateDirectionPlots(app); % Обновляем графики направленности
        end

        % Callback для слайдера угла theta
        function ThetaSliderValueChanged(app, ~)
            powers = [app.PowerSlider1.Value, app.PowerSlider2.Value, app.PowerSlider3.Value];
            frequency = app.FrequencySlider.Value;
            theta = app.ThetaSlider.Value;
            app.appendMessage(['']);
            app.ThetaLabel.Text = ['Угол сканирования = ' num2str(theta, '%.2f') ' градусов'];
            app.Matrix = initializeMatrix(app, powers, frequency, theta);
            updatePlot(app);
            app.appendMessage(['Изменён угол сканирования θ на значение ' num2str(theta, '%.2f') ' градусов']);
            app.appendMessage(['']);
            updateDirectionPlots(app); % Обновляем графики направленности
        end

        % Callback для выпадающего списка ряда Nx
        function RowDropDownValueChanged(app, ~)
            app.SelectedRow = str2double(app.RowDropDown.Value);
            app.appendMessage(['']);
            app.appendMessage(['Выбран ряд Nx = ' num2str(app.SelectedRow)]);
            app.appendMessage(['']);
            updateDirectionPlots(app); % Обновляем графики направленности
        end

        % Callback для выпадающего списка столбца Ny
        function ColumnDropDownValueChanged(app, ~)
            app.SelectedColumn = str2double(app.ColumnDropDown.Value);
            app.appendMessage(['']);
            app.appendMessage(['Выбран столбец Ny = ' num2str(app.SelectedColumn)]);
            updateDirectionPlots(app); % Обновляем графики направленности
        end

        % Callback для кнопки отладки ЭМП
        function DirectionButtonPushed(app, ~)
            app.appendMessage('Запрос параметров ЭМП');
            frequency = app.FrequencySlider.Value;
            powers = [app.PowerSlider1.Value, app.PowerSlider2.Value, app.PowerSlider3.Value];
            
            % Вычисляем параметры для отладочного вывода
            c = 3e8; % Скорость света, м/с
            f = frequency * 1e6; % Частота в Гц
            wavelength = c / f; % Длина волны, м
            k = 2 * pi / wavelength; % Волновое число, радиан/м
            
            % Максимум интенсивности в дальней зоне (используем Fn_3d)
            [~, ~, ~, ~, ~, ~, ~, ~, Fn_3d, ~, ~, ~] = computeAmplitudePattern(app, frequency);
            max_intensity_dB = 10 * log10(max(Fn_3d(:))); % Максимум в дБ
            
            % Выводим информацию
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
            app.appendMessage(['Изменён режим на ' app.ModeSwitch.Value]);
            app.Matrix = initializeMatrix(app, powers, frequency, theta);
            updatePlot(app);
            updateDirectionPlots(app); % Обновляем графики направленности
        end
        
        % Callback для кнопки построения диаграммы направленности
        function DNPButtonPushed(app, ~)
            app.appendMessage('Построение характеристик направленности...');
            updateDirectionPlots(app); % Обновляем графики направленности
        end
    end

    methods (Access = public)

        % Конструктор приложения
        function app = DipoleMatrixApp12
            % Создаем окно
            app.UIFigure = uifigure('Name', 'Модель нагревного стенда СУРА');
            app.UIFigure.Position = [100 100 1366 865];
        
            % Создаём сетку 2×3 для разделения интерфейса
            app.GridLayout = uigridlayout(app.UIFigure, [2 3]);
            app.GridLayout.RowHeight = {'1x', '1x'};
            app.GridLayout.ColumnWidth = {'1x', '3x', '5x'};
            app.GridLayout.RowSpacing = 5; % Расстояние между строками
            app.GridLayout.ColumnSpacing = 5; % Расстояние между столбцами
            app.GridLayout.BackgroundColor = [0.7 0.7 0.7]; % Светло-серая линия разделения
        
            % Левая верхняя часть: Слайдеры и кнопки
            controlsPanel = uipanel(app.GridLayout);
            controlsPanel.Layout.Row = 1;
            controlsPanel.Layout.Column = 1;
            controlsPanel.Title = 'Блок управления стендом СУРА';
            controlsPanel.BackgroundColor = [0.95 0.95 0.95];
        
            % Слайдеры и метки для секции 1
            app.PowerLabel1 = uilabel(controlsPanel);
            app.PowerLabel1.Position = [20 360 200 22];
            app.PowerLabel1.Text = 'Мощность ПКВ-250-1, кВт';
            
            app.PowerSlider1 = uislider(controlsPanel);
            app.PowerSlider1.Limits = [0 300];
            app.PowerSlider1.Value = 0;
            app.PowerSlider1.Position = [20 350 200 3];
            app.PowerSlider1.MajorTicks = [0:50:300];
            app.PowerSlider1.ValueChangedFcn = @(src, event) app.PowerSlider1ValueChanged;
        
            % Слайдеры и метки для секции 2
            app.PowerLabel2 = uilabel(controlsPanel);
            app.PowerLabel2.Position = [20 295 200 22];
            app.PowerLabel2.Text = 'Мощность ПКВ-250-2, кВт';
            
            app.PowerSlider2 = uislider(controlsPanel);
            app.PowerSlider2.Limits = [0 300];
            app.PowerSlider2.Value = 0;
            app.PowerSlider2.Position = [20 285 200 3];
            app.PowerSlider2.MajorTicks = [0:50:300];
            app.PowerSlider2.ValueChangedFcn = @(src, event) app.PowerSlider2ValueChanged;
        
            % Слайдеры и метки для секции 3
            app.PowerLabel3 = uilabel(controlsPanel);
            app.PowerLabel3.Position = [20 230 200 30];
            app.PowerLabel3.Text = 'Мощность ПКВ-250-3, кВт';
            
            app.PowerSlider3 = uislider(controlsPanel);
            app.PowerSlider3.Limits = [0 300];
            app.PowerSlider3.Value = 0;
            app.PowerSlider3.Position = [20 220 200 22];
            app.PowerSlider3.MajorTicks = [0:50:300];
            app.PowerSlider3.ValueChangedFcn = @(src, event) app.PowerSlider3ValueChanged;
        
            % Единый слайдер и метка для частоты
            app.FrequencyLabel = uilabel(controlsPanel);
            app.FrequencyLabel.Position = [20 170 200 22];
            app.FrequencyLabel.Text = 'Частота передатчиков, МГц';
            
            app.FrequencySlider = uislider(controlsPanel);
            app.FrequencySlider.Limits = [3 10];
            app.FrequencySlider.Value = 5;
            app.FrequencySlider.Position = [20 160 200 3];
            app.FrequencySlider.MajorTicks = [3:1:10];
            app.FrequencySlider.ValueChangedFcn = @(src, event) app.FrequencySliderValueChanged;
            
            % Слайдер и метка для угла theta
            app.ThetaLabel = uilabel(controlsPanel);
            app.ThetaLabel.Position = [20 110 200 22];
            app.ThetaLabel.Text = 'Угол луча, градусы';
            
            app.ThetaSlider = uislider(controlsPanel);
            app.ThetaSlider.Limits = [-40 40];
            app.ThetaSlider.Value = 0; % Начальное значение 0 градусов
            app.ThetaSlider.Position = [20 100 200 3];
            app.ThetaSlider.MajorTicks = [-40:10:40];
            app.ThetaSlider.ValueChangedFcn = @(src, event) app.ThetaSliderValueChanged;
            
            % Выпадающий список для выбора ряда Nx
            app.RowLabel = uilabel(controlsPanel);
            app.RowLabel.Position = [8 45 200 22];
            app.RowLabel.Text = 'Выбор ряда Nx';
            
            app.RowDropDown = uidropdown(controlsPanel);
            app.RowDropDown.Items = string(1:12); % Список рядов от 1 до 12
            app.RowDropDown.Value = '1'; % Начальное значение
            app.RowDropDown.Position = [8 25 90 22];
            app.RowDropDown.ValueChangedFcn = @(src, event) app.RowDropDownValueChanged;
            
            % Выпадающий список для выбора столбца Ny
            app.ColumnLabel = uilabel(controlsPanel);
            app.ColumnLabel.Position = [128 45 200 22];
            app.ColumnLabel.Text = 'Выбор столбца Ny';
            
            app.ColumnDropDown = uidropdown(controlsPanel);
            app.ColumnDropDown.Items = string(1:12); % Список столбцов от 1 до 12
            app.ColumnDropDown.Value = '1'; % Начальное значение
            app.ColumnDropDown.Position = [128 25 90 22];
            app.ColumnDropDown.ValueChangedFcn = @(src, event) app.ColumnDropDownValueChanged;
            
            % Тумблер для выбора режима (X-мода или O-мода)
            app.ModeSwitch = uiswitch(controlsPanel, 'slider');
            app.ModeSwitch.Position = [88 -5 90 22];
            app.ModeSwitch.Items = {'X-мода', 'O-мода'};
            app.ModeSwitch.Value = 'X-мода'; % Начальное значение
            app.ModeSwitch.ValueChangedFcn = @(src, event) app.ModeSwitchValueChanged;
        
            % Кнопка для отладки ЭМП
            app.DirectionButton = uicontrol(controlsPanel, 'Style', 'pushbutton');
            app.DirectionButton.String = 'Отладочный вывод ЭМП';
            app.DirectionButton.Position = [20 -30 90 22];
            app.DirectionButton.Callback = @(src, event) app.DirectionButtonPushed(app);
            
            % Кнопка для построения диаграммы направленности
            app.DNPButton = uicontrol(controlsPanel, 'Style', 'pushbutton');
            app.DNPButton.String = 'Построение ДН';
            app.DNPButton.Position = [130 -30 90 22];
            app.DNPButton.Callback = @(src, event) app.DNPButtonPushed(app);
        
            % Средняя верхняя часть: График "Амплитуда диполей по секциям"
            intensityPanel = uipanel(app.GridLayout);
            intensityPanel.Layout.Row = 1;
            intensityPanel.Layout.Column = 2;
            intensityPanel.Title = 'Блок управления и мониторинга параметров сигналов излучателей';
            intensityPanel.BackgroundColor = [0.95 0.95 0.95];
            
            app.UIAxes = uiaxes(intensityPanel);
            app.UIAxes.Position = [10 40 270 270]; % Делаем квадратным: ширина = высота = 270
            
            % Добавляем кнопку "Ввод данных" под тепловой картой
            dataButton = uibutton(intensityPanel, 'Text', 'Ввод данных');
            dataButton.Position = [25 140 100 22];
            dataButton.ButtonPushedFcn = @(src, event) app.OpenMatrixEditor();
            
            % Добавляем цветовую шкалу как отдельный график в intensityPanel
            app.ScaleAxes = uiaxes(intensityPanel);
            app.ScaleAxes.Position = [290 40 40 270]; % Позиция справа от UIAxes: [x y ширина высота]
            colormap(app.ScaleAxes, 'jet');
            % Создаём градиент от 0 до 1
            gradient = linspace(0, 1, 270)'; % 270 значений от 0 до 1 (соответствует высоте шкалы)
            image(app.ScaleAxes, [1 2], [0 1], gradient, 'CDataMapping', 'scaled'); % Use [0 1] for y-axis to match gradient range
            set(app.ScaleAxes, 'YDir', 'normal', 'XTick', [], 'YTick', [0 0.25 0.5 0.75 1], ...
                               'YTickLabel', {'0', '~0.885', '~1.77', '~2.655', '3.54'}, ...
                               'TickLabelInterpreter', 'latex', 'CLim', [0 1]);
            axis(app.ScaleAxes, 'off'); % Убираем видимые оси
            
            % Добавляем контур вокруг шкалы
            % Match the rectangle to the axes position in data coordinates
            rectX = [0 0 3 3 0.5]; % X coordinates for rectangle (adjusted to fit image)
            rectY = [0 1 1 0 0]; % Y coordinates for rectangle (0 to 1 to match gradient)
            hold(app.ScaleAxes, 'on');
            plot(app.ScaleAxes, rectX, rectY, 'k-', 'LineWidth', 1.5);
            hold(app.ScaleAxes, 'off');
            
            % Добавляем пять рисок справа от шкалы
            tickLength = 2; % Length of ticks in data units (relative to [0 1] scale)
            yTickPositions = [0 0.25 0.5 0.75 1]; % Normalized positions
            for i = 1:5
                line(app.ScaleAxes, [2 2 + tickLength], [yTickPositions(i) yTickPositions(i)], ...
                     'Color', 'k', 'LineWidth', 1.5); % Riscs extend to the right
            end
            
            % Добавляем надписи напротив рисок
            labels = {'0', '1.7', '2.6', '3.1', '3.6'};
            for i = 1:5
                text(app.ScaleAxes, 2 + tickLength + 0.05, yTickPositions(i), labels{i}, ...
                     'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
                     'FontSize', 14, 'Interpreter', 'latex');
            end
            
            % Добавляем подпись оси
            ylabel(app.ScaleAxes, 'Мощность (кВт)', 'Rotation', 0, 'HorizontalAlignment', 'right', ...
                   'Position', [-0.1, 0.5, 0]); % Center the label vertically
            % Добавляем надпись "Палитра мощности на эффективных параметрах, кВт"
            title(app.ScaleAxes, 'Палитра мощности на эффективных параметрах, кВт', ...
                  'FontSize', 14, 'Interpreter', 'latex', 'Position', [7, 0.5, 0], 'Rotation', 90);


            % Левая нижняя часть: Центр контроля и отладки
            debugPanel = uipanel(app.GridLayout);
            debugPanel.Layout.Row = 2;
            debugPanel.Layout.Column = [1 2];
            debugPanel.Title = 'Центр контроля и отладки';
            debugPanel.BackgroundColor = [0.95 0.95 0.95];
            
            app.TextArea = uitextarea(debugPanel);
            app.TextArea.Position = [10 10 580 150];
            app.TextArea.Editable = 'off';
        
            % Правая часть: Графики направленности
            directionPanel = uipanel(app.GridLayout);
            directionPanel.Layout.Row = [1 2];
            directionPanel.Layout.Column = 3;
            directionPanel.Title = 'Поле в дальней зоне';
            directionPanel.BackgroundColor = [0.95 0.95 0.95];
            
            % Создаём внутреннюю сетку для графиков направленности
            directionGrid = uigridlayout(directionPanel, [2 3]);
            directionGrid.RowHeight = {'1x', '1x'};
            directionGrid.ColumnWidth = {'1x', '1x', '1x'};
            directionGrid.RowSpacing = 5;
            directionGrid.ColumnSpacing = 5;
            
            % Создаём оси для графиков направленности
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
        
            % Инициализация
            startupFcn(app);
        end
    end

    % Запуск приложения
    methods (Static)
        function run()
            app = DipoleMatrixApp12;
            disp('Интерфейс запущен.');
        end
    end
end