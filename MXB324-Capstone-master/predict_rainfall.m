function rainfall = predict_rainfall(years, state, d_year, f_year)
%PREDICT_RAINFALL Predicts the rainfall per day for the given scenario and
%number of years.
%
%Outputs: 
%- rainfall = The amount of rainfall (mm) that we expect would fall on the
%             each day after t=0. This is a function of days.
%
%Inputs:
%- years = How many years that must be evaluated.
%- state = String input used to specify the climate behaviour being
%          predicted. 
%            - 'drought': A simulation only including drought periods
%            - 'flood': A simulation only including flood periods
%            - 'combination': A simulation only including both drought and flood 
%                             periods
%- d_year = Matrix containing the starting and ending years of drought periods 
%           (1st column = start year, 2nd column = end year).
%- f_year = Matrix containing the starting and ending years of flood periods.
%           (1st column = start year, 2nd column = end year).
%
%Task name: predict_rainfall
%Authors: Caitlin, Jasmine
%
%Revisions [Date, people, reason]:
%- 25/9/2019, Caitlin, implementing drought/flood functionality
%- 29/9/2019, Caitlin, updated so that patterns containing only
%                      drought/flood periods are possible 

%% Handle the input parameters appropriately

% If only the number of years is provided, predict normal weather patterns
% (i.e. no drought/flood periods)
if nargin == 1
    d_year = [];
    f_year = [];
end

% If only flood or only drought behaviour is desired (only one set of year 
% periods is provided). The state input specifies if the input year periods 
% are to be flood or drought. 
if nargin == 3 && strcmp(state,'drought')
    f_year = [];
else if nargin == 3 && strcmp(state,'flood')
        f_year = d_year;
        d_year = [];
    end
end
%% Calculate rf :
rf = 7;
drought_factor = 5/7;
flood_factor = 10/7;

rainfall_function = @(rf, t) rf + rf .* cos(2*pi*t / 365);

%% Predict the rainfall   
rainfall = rainfall_function(rf, 1:(years * 365));

% Implement the drought periods by decreasing the rainfall distribution
for i = 1:length(d_year)
    drought_period = (365*d_year(i, 1)):(365*d_year(i, 2));
    rainfall(drought_period) = rainfall(drought_period) * drought_factor;
end
% Implement the flood periods by increasing the rainfall distribution
for j = 1:length(f_year)
    flood_period = (365*f_year(j, 1)):(365*f_year(j, 2));
    rainfall(flood_period) = rainfall(flood_period) * flood_factor;
end

%% Plot the rainfall

figure();
plot(1:(365 * years), rainfall); 
xlabel('Days');
ylabel('Millimetres of rainfall');

end