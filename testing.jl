using Parameters

mutable struct Person
    name::String
    age::Int 

    function Person(name::String, age::Int = 0)
        new(name, age)
    end
end