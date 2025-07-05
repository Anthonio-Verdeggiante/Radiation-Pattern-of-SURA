classdef DipoleMatrixApp6 < matlab.apps.AppBase

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
            dipole_length = 34; % Длина диполя, м (изменено с 25 м на 34 м)
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
            % Секция 1: столбцы 1-4
            amplitude_variation = 1 + 0.1 * rand(12, 4);
            phase_base = 2 * pi * frequency * (dipole_length / c) * rand(12, 4); % Базовая фаза для первого диполя
            if strcmp(app.ModeSwitch.Value, 'X-мода')
                phase_1 = phase_base; % Фаза первого диполя
                phase_2 = phase_base; % Фаза второго диполя (та же, что у первого)
            else % O-мода
                phase_1 = phase_base; % Фаза первого диполя
                phase_2 = phase_base + pi/2; % Фаза второго диполя (сдвиг на pi/2 вместо pi/4)
            end
            base_matrix(:, 1:4, 1) = sqrt(powers(1)) * amplitude_factor * amplitude_variation .* exp(1i * phase_1);
            base_matrix(:, 1:4, 2) = sqrt(powers(1)) * amplitude_factor * amplitude_variation .* exp(1i * phase_2);
            
            % Секция 2: столбцы 5-8
            amplitude_variation = 1 + 0.1 * rand(12, 4);
            phase_base = 2 * pi * frequency * (dipole_length / c) * rand(12, 4);
            if strcmp(app.ModeSwitch.Value, 'X-мода')
                phase_1 = phase_base;
                phase_2 = phase_base;
            else % O-мода
                phase_1 = phase_base;
                phase_2 = phase_base + pi/2; % Фаза второго диполя (сдвиг на pi/2 вместо pi/4)
            end
            base_matrix(:, 5:8, 1) = sqrt(powers(2)) * amplitude_factor * amplitude_variation .* exp(1i * phase_1);
            base_matrix(:, 5:8, 2) = sqrt(powers(2)) * amplitude_factor * amplitude_variation .* exp(1i * phase_2);
            
            % Секция 3: столбцы 9-12
            amplitude_variation = 1 + 0.1 * rand(12, 4);
            phase_base = 2 * pi * frequency * (dipole_length / c) * rand(12, 4);
            if strcmp(app.ModeSwitch.Value, 'X-мода')
                phase_1 = phase_base;
                phase_2 = phase_base;
            else % O-мода
                phase_1 = phase_base;
                phase_2 = phase_base + pi/2; % Фаза второго диполя (сдвиг на pi/2 вместо pi/4)
            end
            base_matrix(:, 9:12, 1) = sqrt(powers(3)) * amplitude_factor * amplitude_variation .* exp(1i * phase_1);
            base_matrix(:, 9:12, 2) = sqrt(powers(3)) * amplitude_factor * amplitude_variation .* exp(1i * phase_2);
            
            % Учёт взаимодействия (упрощённая модель)
            for i = 1:12
                for j = 1:12
                    for dipole_idx = 1:2 % Для каждого диполя в элементе
                        if abs(base_matrix(i, j, dipole_idx)) > 0
                            coupling_effect = 0;
                            for m = 1:12
                                for n = 1:12
                                    for other_dipole_idx = 1:2
                                        if (m ~= i || n ~= j || dipole_idx ~= other_dipole_idx) && abs(base_matrix(m, n, other_dipole_idx)) > 0
                                            % Расстояние между диполями (i,j) и (m,n)
                                            dist = sqrt(((i-m)*dx_meters)^2 + ((j-n)*dx_meters)^2);
                                            if dist > 0
                                                % Упрощённая модель взаимодействия: затухание ~ 1/r, фазовый сдвиг ~ k*r
                                                coupling = (1/dist) * exp(1i * (2 * pi / wavelength) * dist);
                                                coupling_effect = coupling_effect + coupling * base_matrix(m, n, other_dipole_idx);
                                            end
                                        end
                                    end
                                end
                            end
                            % Корректируем амплитуду и фазу
                            matrix(i, j, dipole_idx) = base_matrix(i, j, dipole_idx) + 0.1 * coupling_effect;
                        end
                    end
                end
            end
        end

        % Функция для вычисления амплитудной характеристики направленности
        function [psi_base, Fn_psi, theta_x, Fn_theta_x, theta_y, Fn_theta_y, theta_x_3d, theta_y_3d, Fn_3d, x_3d, y_3d, Fn_xy] = computeAmplitudePattern(app, frequency)
            % Параметры
            Nx_total = 12; % Число вибраторов вдоль оси x (всего)
            Ny_total = 12; % Число вибраторов вдоль оси y (всего)
            Nx = app.SelectedRow; % Выбранный ряд Nx (индекс)
            Ny = app.SelectedColumn; % Выбранный столбец Ny (индекс)
            d = 25; % Расстояние между излучателями, м
            c = 3e8; % Скорость света, м/с
            ksi = 1; % Коэффициент замедления волны (в радианах, как угол)
            theta_slider = deg2rad(app.ThetaSlider.Value); % Текущее значение theta (в радианах) из слайдера
            theta_offset = 34; % Поправка в 34 градуса 
            
            % Вычисление волнового числа
            f = frequency * 1e6; % Частота в Гц (перевод из МГц)
            wavelength = c / f; % Длина волны, м
            k = 2 * pi / wavelength; % Волновое число, радиан/м
            
            % 1. Зависимость Fn от psi (обобщённая угловая переменная)
            psi_0 = (Nx_total * k * d / 2) * cos(theta_slider + theta_offset - ksi); % Смещение psi
            psi_base = linspace(-10*pi, 10*pi, 1000); % Базовый диапазон psi (фиксированный для графика)
            psi_shifted = psi_base + psi_0; % Смещённый psi для вычисления Fn
            Fn_psi = abs(sin(psi_shifted) ./ (Nx_total * sin(psi_shifted / Nx_total)));
            Fn_psi(isnan(Fn_psi) | isinf(Fn_psi)) = 1; % Устраняем неопределённости
            
            % 2. Зависимость Fn от theta_x и theta_y (для 3D-графика)
            theta_x_3d = linspace(-pi/2, pi/2, 361); % Углы от -90 до 90 градусов (в радианах)
            theta_y_3d = linspace(-pi/2, pi/2, 361); % Углы от -90 до 90 градусов (в радианах)
            [theta_x_3d, theta_y_3d] = meshgrid(theta_x_3d, theta_y_3d);
            Fn_3d = zeros(size(theta_x_3d)); % Инициализируем массив для Fn
            
            % Вычисляем координаты элементов решётки
            x_positions = d * ((1:Nx_total) - (Nx_total + 1)/2); % Координаты x (центрированы вокруг 0)
            y_positions = d * ((1:Ny_total) - (Ny_total + 1)/2); % Координаты y (центрированы вокруг 0)
            
            % Суммируем вклады всех элементов (с учётом двух диполей)
            for i = 1:Nx_total
                for j = 1:Ny_total
                    for dipole_idx = 1:2
                        A_ij = abs(app.Matrix(i, j, dipole_idx)); % Амплитуда элемента (i,j) для диполя
                        phi_ij = angle(app.Matrix(i, j, dipole_idx)); % Фаза элемента (i,j) для диполя
                        % Угол ориентации диполя (45° для первого диполя, 135° для второго)
                        dipole_angle = (dipole_idx == 1) * pi/4 + (dipole_idx == 2) * (pi/4 + pi/2);
                        % Фазовый сдвиг с учётом ориентации диполя
                        phase = k * (x_positions(i) * cos(theta_x_3d + theta_offset - theta_slider - ksi) + ...
                                    y_positions(j) * cos(theta_y_3d + theta_offset - ksi));
                        % Учитываем ориентацию диполя (проекция на направление theta_x, theta_y)
                        directional_factor = cos(dipole_angle - atan2(cos(theta_y_3d), cos(theta_x_3d)));
                        Fn_3d = Fn_3d + A_ij * directional_factor .* exp(1i * (phi_ij + phase));
                    end
                end
            end
            Fn_3d = abs(Fn_3d); % Модуль комплексной суммы
            max_Fn_3d = max(Fn_3d(:));
            if max_Fn_3d > 0
                Fn_3d = Fn_3d / max_Fn_3d; % Нормируем на максимум
            else
                Fn_3d = zeros(size(Fn_3d)); % Если все значения нулевые, устанавливаем нули
            end
            
            % 2. Зависимость Fn от theta_x (срез из Fn_3d для выбранного столбца Ny)
            theta_x = linspace(-pi/2, pi/2, 181); % Углы от -90 до 90 градусов (в радианах)

            % Линейно отображаем Ny (1:12) на индексы theta_y_3d (1:50)
            ny_idx = round(1 + (Ny - 1) * (50 - 1) / (12 - 1)); % Индекс для выбранного столбца Ny
            ny_idx = max(1, min(50, ny_idx)); % Убедимся, что индекс в пределах допустимого

            % Берем срез из Fn_3d вдоль фиксированного theta_y (по столбцу Ny)
            Fn_theta_x = interp1(theta_x_3d(1, :), Fn_3d(ny_idx, :), theta_x, 'linear'); % Интерполяция для более точного соответствия
            Fn_theta_x(isnan(Fn_theta_x)) = 0; % Заменяем NaN на 0
            max_Fn_theta_x = max(Fn_theta_x(:));
            if max_Fn_theta_x > 0
                Fn_theta_x = Fn_theta_x / max_Fn_theta_x; % Нормируем на максимум
            else
                Fn_theta_x = zeros(size(Fn_theta_x)); % Если все значения нулевые, устанавливаем нули
            end

            % 3. Зависимость Fn от theta_y (срез из Fn_3d для выбранного ряда Nx)
            theta_y = linspace(-pi/2, pi/2, 181); % Углы от -90 до 90 градусов (в радианах)

            % Линейно отображаем Nx (1:12) на индексы theta_x_3d (1:50)
            nx_idx = round(1 + (Nx - 1) * (50 - 1) / (12 - 1)); % Индекс для выбранного ряда Nx
            nx_idx = max(1, min(50, nx_idx)); % Убедимся, что индекс в пределах допустимого

            % Берем срез из Fn_3d вдоль фиксированного theta_x (по ряду Nx)
            Fn_theta_y = interp1(theta_y_3d(:, 1), Fn_3d(:, nx_idx), theta_y, 'linear'); % Интерполяция для более точного соответствия
            Fn_theta_y(isnan(Fn_theta_y)) = 0; % Заменяем NaN на 0
            max_Fn_theta_y = max(Fn_theta_y(:));
            if max_Fn_theta_y > 0
                Fn_theta_y = Fn_theta_y / max_Fn_theta_y; % Нормируем на максимум
            else
                Fn_theta_y = zeros(size(Fn_theta_y)); % Если все значения нулевые, устанавливаем нули
            end
            
            % 5. Трёхмерное отображение Fn(x, y) в декартовых координатах
            % Используем тороидальную систему координат
            R = 1; % Радиус большого круга тора
            r = 0.5; % Радиус малого круга тора
            x_3d = (R + r * cos(theta_y_3d)) .* cos(theta_x_3d);
            y_3d = (R + r * cos(theta_y_3d)) .* sin(theta_x_3d);
            
            % Учитываем Якобиан перехода
            jacobian = 1 ./ (sin(theta_x_3d).^2 .* sin(theta_y_3d).^2);
            jacobian(abs(jacobian) < eps) = eps; % Избегаем деления на ноль
            Fn_xy = Fn_3d .* jacobian; % Умножаем на Якобиан
            Fn_xy(isnan(Fn_xy)) = 0; % Заменяем NaN на 0
            max_Fn_xy = max(Fn_xy(:));
            if max_Fn_xy > 0
                Fn_xy = Fn_xy / max_Fn_xy; % Нормируем на максимум
            else
                Fn_xy = zeros(size(Fn_xy)); % Если все значения нулевые, устанавливаем нули
            end

            % Ограничиваем диапазон x, y для читаемости графика
            mask = (abs(x_3d) <= 5) & (abs(y_3d) <= 5);
            x_3d(~mask) = NaN;
            y_3d(~mask) = NaN;
            Fn_xy(~mask) = NaN;

            % Отладочный вывод
            app.appendMessage(['Частота: ' num2str(frequency) ' МГц']);
            app.appendMessage(['Длина волны: ' num2str(wavelength) ' м']);
            app.appendMessage(['Волновое число k: ' num2str(k) ' радиан/м']);
            app.appendMessage(['Смещение psi_0: ' num2str(psi_0) ' радиан']);
            app.appendMessage(['Максимальное значение Fn(psi): ' num2str(max(Fn_psi))]);
            app.appendMessage(['Максимальное значение Fn(theta_x) после нормировки: ' num2str(max_Fn_theta_x)]);
            app.appendMessage(['Максимальное значение Fn(theta_y) после нормировки: ' num2str(max_Fn_theta_y)]);
            app.appendMessage(['Максимальное значение Fn_3d после нормировки: ' num2str(max_Fn_3d)]);
        end

        % Обновление графика "Амплитуда диполей по секциям"
        function updatePlot(app)
            % Вычисляем среднюю амплитуду по обоим диполям для отображения
            amplitude = mean(abs(app.Matrix), 3); % Средняя амплитуда по третьему измерению (два диполя)
            if isempty(app.UIAxes) || ~isvalid(app.UIAxes)
                app.UIAxes = uiaxes(app.GridLayout);
                app.UIAxes.Layout.Row = 1;
                app.UIAxes.Layout.Column = 2;
            end
            imagesc(app.UIAxes, amplitude);
            colormap(app.UIAxes, 'jet');
            colorbar(app.UIAxes);
            axis(app.UIAxes, [0.5 11.5 0.5 11.5]);
            axis(app.UIAxes, 'equal');
            
            % Добавляем сетку
            hold(app.UIAxes, 'on');
            for k = 0.5:1:13
                plot(app.UIAxes, [k k], [0 12.5], 'k-', 'LineWidth', 0.5);
                plot(app.UIAxes, [0 12.5], [k k], 'k-', 'LineWidth', 0.5);
            end
            hold(app.UIAxes, 'off');
            title(app.UIAxes, 'Амплитуда диполей по секциям');
            drawnow;
        end

        % Обновление графиков направленности
        function updateDirectionPlots(app)
            frequency = app.FrequencySlider.Value;
            [psi_base, Fn_psi, theta_x, Fn_theta_x, theta_y, Fn_theta_y, theta_x_3d, theta_y_3d, Fn_3d, x_3d, y_3d, Fn_xy] = computeAmplitudePattern(app, frequency);
            
            % График 1: Fn(psi)
            cla(app.AxesPsi); % Очищаем оси перед перерисовкой
            plot(app.AxesPsi, psi_base, Fn_psi, 'b-');
            xlabel(app.AxesPsi, '\psi (радианы)');
            ylabel(app.AxesPsi, 'F_n(\psi)');
            title(app.AxesPsi, 'Амплитудная характеристика от \psi');
            grid(app.AxesPsi, 'on');
            
            % График 2: Fn(theta_x) (вдоль оси x)
            cla(app.AxesThetaX); % Очищаем оси перед перерисовкой
            plot(app.AxesThetaX, rad2deg(theta_x), Fn_theta_x, 'r-');
            xlabel(app.AxesThetaX, '\theta_x (градусы)');
            ylabel(app.AxesThetaX, 'F_n(\theta_x) (нормализованная)');
            title(app.AxesThetaX, ['Амплитудная характеристика от \theta_x (ряд ' num2str(app.SelectedRow) ')']);
            grid(app.AxesThetaX, 'on');
            
            % График 3: Fn(theta_y) (вдоль оси y)
            cla(app.AxesThetaY); % Очищаем оси перед перерисовкой
            plot(app.AxesThetaY, rad2deg(theta_y), Fn_theta_y, 'g-');
            xlabel(app.AxesThetaY, '\theta_y (градусы)');
            ylabel(app.AxesThetaY, 'F_n(\theta_y) (нормализованная)');
            title(app.AxesThetaY, ['Амплитудная характеристика от \theta_y (столбец ' num2str(app.SelectedColumn) ')']);
            grid(app.AxesThetaY, 'on');
            
            % График 4: 3D отображение Fn(theta_x, theta_y) для всей решётки
            cla(app.AxesTheta3D); % Очищаем оси перед перерисовкой
            surf(app.AxesTheta3D, rad2deg(theta_x_3d), rad2deg(theta_y_3d), Fn_3d, 'EdgeColor', 'interp');
            xlabel(app.AxesTheta3D, '\theta_x (градусы)');
            ylabel(app.AxesTheta3D, '\theta_y (градусы)');
            zlabel(app.AxesTheta3D, 'F_n(\theta_x, \theta_y)');
            title(app.AxesTheta3D, '3D характеристика F_n(\theta_x, \theta_y) (вся решётка)');
            colormap(app.AxesTheta3D, 'jet');
            colorbar(app.AxesTheta3D);
            grid(app.AxesTheta3D, 'on');
            
            % График 5: 3D отображение Fn(x, y)
            cla(app.AxesXY); % Очищаем оси перед перерисовкой
            surf(app.AxesXY, x_3d, y_3d, Fn_xy, 'EdgeColor', 'interp');
            xlabel(app.AxesXY, 'x');
            ylabel(app.AxesXY, 'y');
            zlabel(app.AxesXY, 'F_n(x, y)');
            title(app.AxesXY, '3D характеристика F_n(x, y) (вся решётка)');
            colormap(app.AxesXY, 'jet');
            colorbar(app.AxesXY);
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
            app.appendMessage(['Изменена частота на ' num2str(frequency, '%.2f') ' МГц']);
            updateDirectionPlots(app); % Обновляем графики направленности
        end

        % Callback для слайдера угла theta
        function ThetaSliderValueChanged(app, ~)
            powers = [app.PowerSlider1.Value, app.PowerSlider2.Value, app.PowerSlider3.Value];
            frequency = app.FrequencySlider.Value;
            theta = app.ThetaSlider.Value;
            app.ThetaLabel.Text = ['Угол \theta = ' num2str(theta, '%.2f') ' градусов'];
            app.Matrix = initializeMatrix(app, powers, frequency, theta);
            updatePlot(app);
            app.appendMessage(['Изменён угол \theta на ' num2str(theta, '%.2f') ' градусов']);
            updateDirectionPlots(app); % Обновляем графики направленности
        end

        % Callback для выпадающего списка ряда Nx
        function RowDropDownValueChanged(app, ~)
            app.SelectedRow = str2double(app.RowDropDown.Value);
            app.appendMessage(['Выбран ряд Nx = ' num2str(app.SelectedRow)]);
            updateDirectionPlots(app); % Обновляем графики направленности
        end

        % Callback для выпадающего списка столбца Ny
        function ColumnDropDownValueChanged(app, ~)
            app.SelectedColumn = str2double(app.ColumnDropDown.Value);
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
            app.appendMessage(['Частота излучения: ' num2str(frequency) ' МГц']);
            app.appendMessage(['Длина волны: ' num2str(wavelength) ' м']);
            app.appendMessage(['Волновое число k: ' num2str(k) ' радиан/м']);
            app.appendMessage(['Максимум интенсивности в дальней зоне: ' num2str(max_intensity_dB) ' дБ']);
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
        function app = DipoleMatrixApp6
            % Создаем окно
            app.UIFigure = uifigure('Name', 'Модель нагревного стенда СУРА');
            app.UIFigure.Position = [100 100 1366 865];

            % Создаём сетку 2×3 для разделения интерфейса
            app.GridLayout = uigridlayout(app.UIFigure, [2 3]);
            app.GridLayout.RowHeight = {'1x', '1x'};
            app.GridLayout.ColumnWidth = {'1x', '2x', '5x'};
            app.GridLayout.RowSpacing = 5; % Расстояние между строками
            app.GridLayout.ColumnSpacing = 5; % Расстояние между столбцами
            app.GridLayout.BackgroundColor = [0.7 0.7 0.7]; % Светло-серая линия разделения

            % Левая верхняя часть: Слайдеры и кнопки
            controlsPanel = uipanel(app.GridLayout);
            controlsPanel.Layout.Row = 1;
            controlsPanel.Layout.Column = 1;
            controlsPanel.Title = 'Управление стендом СУРА';
            controlsPanel.BackgroundColor = [0.95 0.95 0.95];

            % Слайдеры и метки для секции 1
            app.PowerLabel1 = uilabel(controlsPanel);
            app.PowerLabel1.Position = [20 360 200 22];
            app.PowerLabel1.Text = 'Мощность ПКВ-250-1, кВт';
            
            app.PowerSlider1 = uislider(controlsPanel);
            app.PowerSlider1.Limits = [0 300];
            app.PowerSlider1.Value = 6.25;
            app.PowerSlider1.Position = [20 350 180 3];
            app.PowerSlider1.MajorTicks = [0:50:300];
            app.PowerSlider1.ValueChangedFcn = @(src, event) app.PowerSlider1ValueChanged;

            % Слайдеры и метки для секции 2
            app.PowerLabel2 = uilabel(controlsPanel);
            app.PowerLabel2.Position = [20 295 200 22];
            app.PowerLabel2.Text = 'Мощность ПКВ-250-2, кВт';
            
            app.PowerSlider2 = uislider(controlsPanel);
            app.PowerSlider2.Limits = [0 300];
            app.PowerSlider2.Value = 6.25;
            app.PowerSlider2.Position = [20 285 180 3];
            app.PowerSlider2.MajorTicks = [0:50:300];
            app.PowerSlider2.ValueChangedFcn = @(src, event) app.PowerSlider2ValueChanged;

            % Слайдеры и метки для секции 3
            app.PowerLabel3 = uilabel(controlsPanel);
            app.PowerLabel3.Position = [20 230 200 30];
            app.PowerLabel3.Text = 'Мощность ПКВ-250-3, кВт';
            
            app.PowerSlider3 = uislider(controlsPanel);
            app.PowerSlider3.Limits = [0 300];
            app.PowerSlider3.Value = 6.25;
            app.PowerSlider3.Position = [20 220 180 22];
            app.PowerSlider3.MajorTicks = [0:50:300];
            app.PowerSlider3.ValueChangedFcn = @(src, event) app.PowerSlider3ValueChanged;

            % Единый слайдер и метка для частоты
            app.FrequencyLabel = uilabel(controlsPanel);
            app.FrequencyLabel.Position = [20 170 200 22];
            app.FrequencyLabel.Text = 'Частота передатчиков, МГц';
            
            app.FrequencySlider = uislider(controlsPanel);
            app.FrequencySlider.Limits = [3 10];
            app.FrequencySlider.Value = 5;
            app.FrequencySlider.Position = [20 160 180 3];
            app.FrequencySlider.MajorTicks = [3:1:10];
            app.FrequencySlider.ValueChangedFcn = @(src, event) app.FrequencySliderValueChanged;
            
            % Слайдер и метка для угла theta
            app.ThetaLabel = uilabel(controlsPanel);
            app.ThetaLabel.Position = [20 110 200 22];
            app.ThetaLabel.Text = 'Угол фазирования, градусы';
            
            app.ThetaSlider = uislider(controlsPanel);
            app.ThetaSlider.Limits = [-40 40];
            app.ThetaSlider.Value = 0; % Начальное значение 0 градусов
            app.ThetaSlider.Position = [20 100 180 3];
            app.ThetaSlider.MajorTicks = [-40:10:40];
            app.ThetaSlider.ValueChangedFcn = @(src, event) app.ThetaSliderValueChanged;
            
            % Выпадающий список для выбора ряда Nx
            app.RowLabel = uilabel(controlsPanel);
            app.RowLabel.Position = [20 45 200 22];
            app.RowLabel.Text = 'Выбор ряда Nx';
            
            app.RowDropDown = uidropdown(controlsPanel);
            app.RowDropDown.Items = string(1:12); % Список рядов от 1 до 12
            app.RowDropDown.Value = '1'; % Начальное значение
            app.RowDropDown.Position = [20 25 90 22];
            app.RowDropDown.ValueChangedFcn = @(src, event) app.RowDropDownValueChanged;
            
            % Выпадающий список для выбора столбца Ny
            app.ColumnLabel = uilabel(controlsPanel);
            app.ColumnLabel.Position = [130 45 200 22];
            app.ColumnLabel.Text = 'Выбор столбца Ny';
            
            app.ColumnDropDown = uidropdown(controlsPanel);
            app.ColumnDropDown.Items = string(1:12); % Список столбцов от 1 до 12
            app.ColumnDropDown.Value = '1'; % Начальное значение
            app.ColumnDropDown.Position = [130 25 90 22];
            app.ColumnDropDown.ValueChangedFcn = @(src, event) app.ColumnDropDownValueChanged;
            
            % Тумблер для выбора режима (X-мода или O-мода)
            app.ModeSwitch = uiswitch(controlsPanel, 'slider');
            app.ModeSwitch.Position = [20 -5 90 22];
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
            intensityPanel.Title = 'Интенсивность излучателей';
            intensityPanel.BackgroundColor = [0.95 0.95 0.95];

            app.UIAxes = uiaxes(intensityPanel);
            app.UIAxes.Position = [10 10 400 300]; % Начальный размер графика

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
            app = DipoleMatrixApp6;
            disp('Интерфейс запущен.');
        end
    end
end