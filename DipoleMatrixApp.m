classdef DipoleMatrixApp < matlab.apps.AppBase

    % Свойства приложения
    properties (Access = public)
        UIFigure          matlab.ui.Figure
        PowerSlider1      matlab.ui.control.Slider % Мощность секции 1
        PowerSlider2      matlab.ui.control.Slider % Мощность секции 2
        PowerSlider3      matlab.ui.control.Slider % Мощность секции 3
        FrequencySlider   matlab.ui.control.Slider % Единый слайдер частоты
        PowerLabel1       matlab.ui.control.Label
        PowerLabel2       matlab.ui.control.Label
        PowerLabel3       matlab.ui.control.Label
        FrequencyLabel    matlab.ui.control.Label
        UIAxes            matlab.ui.control.UIAxes % Для интенсивности матрицы
        Matrix            % Матрица комплексных чисел
        DirectionButton   matlab.ui.control.UIControl % Кнопка для отладки ЭМП
        DNPButton         matlab.ui.control.UIControl % Кнопка для построения ДН
        TextArea          matlab.ui.control.TextArea % Область для сообщений
        DirectionAxes     matlab.ui.Figure % Окно для диаграммы направленности
    end

    methods (Access = private)

        % Инициализация приложения
        function startupFcn(app)
            app.appendMessage('Инициализация приложения...');
            % Изначальная матрица 12x12
            app.Matrix = initializeMatrix(app, [6.25, 6.25, 6.25], app.FrequencySlider.Value);
            updatePlot(app); % Обновляем график
            app.appendMessage('График должен быть инициализирован');
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
        function matrix = initializeMatrix(app, powers, frequency)
            matrix = zeros(12, 12);
            c = 3e8; % Скорость света, м/с
            dipole_length = 25; % Длина диполя, м
            resonant_frequency = c / (2 * dipole_length); % Резонансная частота
            resonant_frequency_MHz = resonant_frequency / 1e6; % Перевод в МГц (~6 МГц)
            
            % Резонансная характеристика (гауссова функция)
            bandwidth = 5.0; % Полоса пропускания, МГц
            amplitude_factor = exp(-((frequency - resonant_frequency_MHz) / bandwidth).^2);
            
            % Шаг между диполями в метрах
            dx = 0.5; % Шаг в долях длины волны
            wavelength = c / (frequency * 1e6);
            dx_meters = dx * wavelength; % Шаг в метрах
            
            % Базовые амплитуды и фазы (без взаимодействия)
            base_matrix = zeros(12, 12, 'like', 1i);
            % Секция 1: столбцы 1-4
            amplitude_variation = 1 + 0.1 * rand(12, 4);
            phase = 2 * pi * frequency * (dipole_length / c) * rand(12, 4);
            base_matrix(:, 1:4) = sqrt(powers(1)) * amplitude_factor * amplitude_variation .* exp(1i * phase);
            
            % Секция 2: столбцы 5-8
            amplitude_variation = 1 + 0.1 * rand(12, 4);
            phase = 2 * pi * frequency * (dipole_length / c) * rand(12, 4);
            base_matrix(:, 5:8) = sqrt(powers(2)) * amplitude_factor * amplitude_variation .* exp(1i * phase);
            
            % Секция 3: столбцы 9-12
            amplitude_variation = 1 + 0.1 * rand(12, 4);
            phase = 2 * pi * frequency * (dipole_length / c) * rand(12, 4);
            base_matrix(:, 9:12) = sqrt(powers(3)) * amplitude_factor * amplitude_variation .* exp(1i * phase);
            
            % Учёт взаимодействия (упрощённая модель)
            for i = 1:12
                for j = 1:12
                    if abs(base_matrix(i, j)) > 0
                        coupling_effect = 0;
                        for m = 1:12
                            for n = 1:12
                                if (m ~= i || n ~= j) && abs(base_matrix(m, n)) > 0
                                    % Расстояние между диполями (i,j) и (m,n)
                                    dist = sqrt(((i-m)*dx_meters)^2 + ((j-n)*dx_meters)^2);
                                    if dist > 0
                                        % Упрощённая модель взаимодействия: затухание ~ 1/r, фазовый сдвиг ~ k*r
                                        coupling = (1/dist) * exp(1i * (2 * pi / wavelength) * dist);
                                        coupling_effect = coupling_effect + coupling * base_matrix(m, n);
                                    end
                                end
                            end
                        end
                        % Корректируем амплитуду и фазу
                        matrix(i, j) = base_matrix(i, j) + 0.1 * coupling_effect;
                    end
                end
            end
        end

        % Новая функция для вычисления амплитудной характеристики направленности
        function [psi, Fn] = computeAmplitudePattern(app, frequency)
            % Параметры
            Nx = 12; % Число вибраторов вдоль оси x
            d = 25; % Расстояние между излучателями, м
            c = 3e8; % Скорость света, м/с
            theta = deg2rad(90); % Направление излучения (90° — строго вверх)
            ksi = 1; % Коэффициент замедления волны (в вакууме)
            
            % Вычисление волнового числа
            f = frequency * 1e6; % Частота в Гц (перевод из МГц)
            wavelength = c / f; % Длина волны, м
            k = 2 * pi / wavelength; % Волновое число, радиан/м
            
            % Диапазон для psi (обобщённая угловая переменная)
            psi = linspace(-10*pi, 10*pi, 1000); % Диапазон от -10π до 10π, 1000 точек
            
            % Вычисление амплитудной характеристики Fn
            Fn = abs(sin(psi) ./ (Nx * sin(psi / Nx)));
            
            % Устранение неопределённостей (где sin(psi/Nx) ≈ 0)
            Fn(isnan(Fn) | isinf(Fn)) = 1; % В точках неопределённости (psi = 0, ±Nxπ, ...) Fn = 1
            
            % Отладочный вывод
            app.appendMessage(['Частота: ' num2str(frequency) ' МГц']);
            app.appendMessage(['Длина волны: ' num2str(wavelength) ' м']);
            app.appendMessage(['Волновое число k: ' num2str(k) ' радиан/м']);
            app.appendMessage(['Максимальное значение Fn: ' num2str(max(Fn))]);
        end

        % Обновление графика
        function updatePlot(app)
            app.appendMessage('Обновление графика...');
            amplitude = abs(app.Matrix); % Амплитуда для отображения
            if isempty(app.UIAxes) || ~isvalid(app.UIAxes)
                app.UIAxes = uiaxes(app.UIFigure, 'Position', [250 180 500 500]);
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

        % Callback для слайдера мощности секции 1
        function PowerSlider1ValueChanged(app, ~)
            powers = [app.PowerSlider1.Value, app.PowerSlider2.Value, app.PowerSlider3.Value];
            frequency = app.FrequencySlider.Value;
            app.PowerLabel1.Text = ['Мощность ПКВ-250-1 = ' num2str(powers(1), '%.2f') ' кВт'];
            app.Matrix = initializeMatrix(app, powers, frequency);
            updatePlot(app);
            app.appendMessage(['Изменена мощность ПКВ-250-1 на ' num2str(powers(1), '%.2f') ' кВт']);
        end

        % Callback для слайдера мощности секции 2
        function PowerSlider2ValueChanged(app, ~)
            powers = [app.PowerSlider1.Value, app.PowerSlider2.Value, app.PowerSlider3.Value];
            frequency = app.FrequencySlider.Value;
            app.PowerLabel2.Text = ['Мощность ПКВ-250-2 = ' num2str(powers(2), '%.2f') ' кВт'];
            app.Matrix = initializeMatrix(app, powers, frequency);
            updatePlot(app);
            app.appendMessage(['Изменена мощность ПКВ-250-2 на ' num2str(powers(2), '%.2f') ' кВт']);
        end

        % Callback для слайдера мощности секции 3
        function PowerSlider3ValueChanged(app, ~)
            powers = [app.PowerSlider1.Value, app.PowerSlider2.Value, app.PowerSlider3.Value];
            frequency = app.FrequencySlider.Value;
            app.PowerLabel3.Text = ['Мощность ПКВ-250-3 = ' num2str(powers(3), '%.2f') ' кВт'];
            app.Matrix = initializeMatrix(app, powers, frequency);
            updatePlot(app);
            app.appendMessage(['Изменена мощность ПКВ-250-3 на ' num2str(powers(3), '%.2f') ' кВт']);
        end

        % Callback для слайдера частоты
        function FrequencySliderValueChanged(app, ~)
            powers = [app.PowerSlider1.Value, app.PowerSlider2.Value, app.PowerSlider3.Value];
            frequency = app.FrequencySlider.Value;
            app.FrequencyLabel.Text = ['Частота передатчиков = ' num2str(frequency, '%.2f') ' МГц'];
            app.Matrix = initializeMatrix(app, powers, frequency);
            updatePlot(app);
            app.appendMessage(['Изменена частота на ' num2str(frequency, '%.2f') ' МГц']);
        end

        % Callback для кнопки отладки ЭМП
        function DirectionButtonPushed(app, ~)
            app.appendMessage('Запрос параметров ЭМП');
            frequency = app.FrequencySlider.Value;
            % Пока просто выводим сообщение, функционал можно доработать позже
            app.appendMessage(['Частота для отладки: ' num2str(frequency) ' МГц']);
            drawnow;
        end

        % Callback для кнопки построения диаграммы направленности
        function DNPButtonPushed(app, ~)
            app.appendMessage('Построение амплитудной характеристики направленности...');
            frequency = app.FrequencySlider.Value;
            [psi, Fn] = computeAmplitudePattern(app, frequency);
            
            % Создание нового интерфейса для диаграммы
            if isempty(app.DirectionAxes) || ~isvalid(app.DirectionAxes)
                app.DirectionAxes = uifigure('Name', 'Амплитудная характеристика направленности', 'Position', [400 100 600 400]);
            end
            
            % Создаём или обновляем оси
            directionAxes = findobj(app.DirectionAxes, 'Type', 'axes');
            if isempty(directionAxes)
                directionAxes = uiaxes(app.DirectionAxes, 'Position', [50 50 500 300]);
            end
            
            % Построение графика
            cla(directionAxes);
            plot(directionAxes, psi, Fn, 'b-');
            xlabel(directionAxes, '\psi (радианы)');
            ylabel(directionAxes, 'F_n(\psi)');
            title(directionAxes, 'Амплитудная характеристика направленности от \psi');
            grid(directionAxes, 'on');
            drawnow;
        end
    end

    methods (Access = public)

        % Конструктор приложения
        function app = DipoleMatrixApp
            % Создаем окно
            app.UIFigure = uifigure('Name', 'Интенсивность излучения в ближней зоне');
            app.UIFigure.Position = [100 100 900 800];

            % Слайдеры и метки для секции 1
            app.PowerLabel1 = uilabel(app.UIFigure);
            app.PowerLabel1.Position = [20 580 200 22];
            app.PowerLabel1.Text = 'Мощность ПКВ-250-1, кВт';
            
            app.PowerSlider1 = uislider(app.UIFigure);
            app.PowerSlider1.Limits = [0 300];
            app.PowerSlider1.Value = 6.25;
            app.PowerSlider1.Position = [20 570 180 3];
            app.PowerSlider1.MajorTicks = [0:50:300];
            app.PowerSlider1.ValueChangedFcn = @(src, event) app.PowerSlider1ValueChanged;

            % Слайдеры и метки для секции 2
            app.PowerLabel2 = uilabel(app.UIFigure);
            app.PowerLabel2.Position = [20 480 200 22];
            app.PowerLabel2.Text = 'Мощность ПКВ-250-2, кВт';
            
            app.PowerSlider2 = uislider(app.UIFigure);
            app.PowerSlider2.Limits = [0 300];
            app.PowerSlider2.Value = 6.25;
            app.PowerSlider2.Position = [20 470 180 3];
            app.PowerSlider2.MajorTicks = [0:50:300];
            app.PowerSlider2.ValueChangedFcn = @(src, event) app.PowerSlider2ValueChanged;

            % Слайдеры и метки для секции 3
            app.PowerLabel3 = uilabel(app.UIFigure);
            app.PowerLabel3.Position = [20 380 200 30];
            app.PowerLabel3.Text = 'Мощность ПКВ-250-3, кВт';
            
            app.PowerSlider3 = uislider(app.UIFigure);
            app.PowerSlider3.Limits = [0 300];
            app.PowerSlider3.Value = 6.25;
            app.PowerSlider3.Position = [20 370 180 22];
            app.PowerSlider3.MajorTicks = [0:50:300];
            app.PowerSlider3.ValueChangedFcn = @(src, event) app.PowerSlider3ValueChanged;

            % Единый слайдер и метка для частоты
            app.FrequencyLabel = uilabel(app.UIFigure);
            app.FrequencyLabel.Position = [20 170 200 22];
            app.FrequencyLabel.Text = 'Частота передатчиков, МГц';
            
            app.FrequencySlider = uislider(app.UIFigure);
            app.FrequencySlider.Limits = [3 10];
            app.FrequencySlider.Value = 5;
            app.FrequencySlider.Position = [20 160 180 3];
            app.FrequencySlider.MajorTicks = [3:0.5:9 9.5 10];
            app.FrequencySlider.ValueChangedFcn = @(src, event) app.FrequencySliderValueChanged;
            
            % Кнопка для отладки ЭМП
            app.DirectionButton = uicontrol(app.UIFigure, 'Style', 'pushbutton');
            app.DirectionButton.String = 'Отладочный вывод ЭМП';
            app.DirectionButton.Position = [20 300 200 30];
            app.DirectionButton.Callback = @(src, event) app.DirectionButtonPushed(app);
            drawnow;

            % Кнопка для построения диаграммы направленности
            app.DNPButton = uicontrol(app.UIFigure, 'Style', 'pushbutton');
            app.DNPButton.String = 'Построение ДН';
            app.DNPButton.Position = [20 230 200 30];
            app.DNPButton.Callback = @(src, event) app.DNPButtonPushed(app);
            drawnow;

            % Область для сообщений
            app.TextArea = uitextarea(app.UIFigure);
            app.TextArea.Position = [20 20 250 70];
            app.TextArea.Editable = 'off';

            % Создаем область для графика
            app.UIAxes = uiaxes(app.UIFigure);
            app.UIAxes.Position = [250 180 500 500];

            % Инициализация
            startupFcn(app);
        end
    end

    % Запуск приложения
    methods (Static)
        function run()
            app = DipoleMatrixApp;
            disp('Интерфейс запущен.');
        end
    end
end