using CSV
import DataFrames as DF
using MLJ

# Load the Titanic dataset
function load_titanic_data(file_path::String)
    df = DF.DataFrame(CSV.File(file_path))

    X = DF.select(df, DF.Not([:Survived, :Name, :Ticket]))
    y = DF.select(df, :Survived)
    return X, y
end

# # Function to preprocess the Titanic dataset
function preprocess_titanic_data(train_X, train_y)
    # println(first(train_X, 10) |> pretty)
    # println(first(train_y, 10) |> pretty)

    train_y2 = train_y[!, 1]  # Convert DataFrame column to vector

    # println(MLJ.models(matching(train_X, train_y2)))
    # println(typeof(train_X))
    # println(typeof(train_y2))

    # println(first(train_X, 10) |> pretty)

    function sex_to_numeric(x)
        if !ismissing(x)
            if x == "female"
                return 1 
            elseif x == "male" 
                return 0 
            end
        else
            return -1  # Assign a numeric value for missing values
        end
    end
    @. train_X.Sex = sex_to_numeric(train_X.Sex)

    function embarked_to_numeric(x)
        if !ismissing(x)
            if x == "S" 
                return 0 
            elseif x == "C" 
                return 1 
            elseif x == "Q" 
                return 2 
            end
        else
            return -1  # Assign a numeric value for missing values
        end
    end
    @. train_X.Embarked = embarked_to_numeric(train_X.Embarked)

    function cabin_to_numeric(x)
        if !ismissing(x)
            return 1
        else
            return 0  # Assign a numeric value for missing values
        end
    end
    @. train_X.Cabin = cabin_to_numeric(train_X.Cabin)

    train_X.Age = coalesce.(train_X.Age, mean(skipmissing(train_X.Age)))  # Fill missing Age with mean
    train_X.Fare = coalesce.(train_X.Fare, mean(skipmissing(train_X.Fare)))  # Fill missing Fare with mean

    # println(first(train_X, 10) |> pretty)

    return train_X, train_y
end

file_path = "train.csv"
train_X, train_y = load_titanic_data(file_path)
train_X, train_y = preprocess_titanic_data(train_X, train_y)

println("Training data loaded and preprocessed successfully.")

println(first(train_X, 10) |> pretty)
# println(schema(train_X) |> pretty)
