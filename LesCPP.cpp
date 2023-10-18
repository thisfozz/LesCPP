#include <iostream>

class Array {
private:
    int* current_array_;
    int size_;

    static int generate_random_number() {
        std::srand(static_cast<unsigned int>(std::time(nullptr)));
        const int random_number = std::rand() % 6 + 5;
        return random_number;
    }

public:
    explicit Array(const int& size) :size_(size)
    {
        current_array_ = new int[size_];
        for (int i = 0; i < size; i++)
        {
            current_array_[i] = generate_random_number();
        }
    }
    Array() = default;

    Array(const Array& other_array) : size_(other_array.size_) {
        current_array_ = new int[size_];
        for (int i = 0; i < size_; i++)
        {
            current_array_[i] = other_array.current_array_[i];
        }
    }

    void ChangeSize(const int& size) {
        const int* tmp_array = current_array_;
        current_array_ = new int[size];

        for (int i = 0; i < size; i++)
        {
            current_array_[i] = tmp_array[i];
        }

        delete[] tmp_array;
        size_ = size;
    }

    void print_array() const
    {
        if (size_ <= 0) return;
        for (size_t i = 0; i < size_; ++i) {
            std::cout << current_array_[i] << " ";
        }
    }
    void sort() const
    {
        if (size_ <= 0) return;
        for (int i = 1; i < size_; i++)
        {
            int j = i - 1;
            while (j >= 0 && current_array_[j] > current_array_[j + 1])
            {
                std::swap(current_array_[j], current_array_[j + 1]);
                j--;
            }
        }
    }

    int find_max() const
    {
        if (size_ <= 0) return -1;
        int max = current_array_[0];
        for (size_t i = 1; i < size_; ++i) {
            if (current_array_[i] > max) {
                max = current_array_[i];
            }
        }
        return max;
    }

    int find_min() const
    {
        if (size_ <= 0) return -1;
        int min = current_array_[0];
        for (size_t i = 1; i < size_; ++i) {
            if (current_array_[i] < min) {
                min = current_array_[i];
            }
        }
        return min;
    }

    int& operator[](const unsigned int index) const
    {
        if (index > size_)
            throw std::exception("Index of Range");
        return current_array_[index];
    }
    int operator()(int value) {
        for (size_t i = 0; i < size_; ++i)
        {
            current_array_[i] + value;
        }
    }

    explicit operator int() const
    {
        int sum = 0;
        for (size_t i = 0; i < size_; ++i)
        {
            sum += current_array_[i];
        }
        return sum;
    }
    operator char* () const
    {
        char* str = new char[size_ + 1];

        for (size_t i = 0; i < size_; ++i)
        {
            str[i] = static_cast<char>(current_array_[i]);
        }
        str[size_] = '\0';

        return str;
    }

    ~Array() {
        delete[] current_array_;
    }
};

template<typename T>
class Var
{
public:
    Var(const T variable) : variable_(variable) {}

    Var operator +(const Var& other)
    {
        Var result;
    	result.variable_ = variable_ + other.variable_;
        return result;
    }
    Var operator +=(const Var& other)
    {
        variable_ += other.variable_;
        return *this;
    }

    Var operator -(const Var& other)
    {
        if (typeid(T) == typeid(std::string) || typeid(T) == typeid(char)) throw std::exception("Invalid argument");
        Var result;
        result.variable_ = variable_ - other.variable_;

        return result;
    }
    Var operator -=(const Var& other)
    {
        if (typeid(T) == typeid(std::string) || typeid(T) == typeid(char)) throw std::exception("Invalid argument");
        variable_ -= other.variable_;
        return *this;
    }

    Var operator *(const Var& other)
    {
        if (typeid(T) == typeid(std::string) || typeid(T) == typeid(char))
        {
            Var result;
            const std::string& str1 = variable_;
            const std::string& str2 = other.variable_;

            for (char c : str1) {
                if (str2.find(c) != std::string::npos) {
                    result += c;
                }
            }

            return result;
        }
        Var result;
        result.variable_ = variable_ * other.variable_;

        return result;
    }
    Var operator *=(const Var& other)
    {
        if (typeid(T) == typeid(std::string) || typeid(T) == typeid(char))
        {
            Var result;
            const std::string& str1 = variable_;
            const std::string& str2 = other.variable_;

            for (char c : str1) {
                if (str2.find(c) != std::string::npos) {
                    result += c;
                }
            }

            return result;
        }
        variable_ *= other.variable_;
        return *this;
    }

    Var operator /(const Var& other)
    {
        if (typeid(T) == typeid(std::string) || typeid(T) == typeid(char))
        {
            Var result;
            const std::string& str1 = variable_;
            const std::string& str2 = other.variable_;

            for (char c : str1) {
                if (str2.find(c) == std::string::npos) {
                    result += c;
                }
            }

            return result;
        }

        if(other.variable_ != 0)
        {
            Var result;
            result.variable_ = variable_ / other.variable_;
            return result;
        }else
        {
            throw std::exception("Divide by zero");
        }
    }
    Var operator /=(const Var& other)
    {
        if (typeid(T) == typeid(std::string) || typeid(T) == typeid(char))
        {
            Var result;
            const std::string& str1 = variable_;
            const std::string& str2 = other.variable_;

            for (char c : str1) {
                if (str2.find(c) == std::string::npos) {
                    result += c;
                }
            }

            return result;
        }

        if (other.variable_ != 0)
        {
            variable_ /= other.variable_;
            return *this;
        }
        else throw std::exception("Divide by zero");
    }

    bool operator <(const Var& other)
    {
	    if(typeid(T) == typeid(std::string) || typeid(T) == typeid(char))
	    {
            return variable_.length() < other.variable_.length();
	    }else
	    {
            return variable_ < other.variable_;
	    }
    }
    bool operator <=(const Var& other)
    {
        if (typeid(T) == typeid(std::string) || typeid(T) == typeid(char))
        {
            return variable_.length() <= other.variable_.length();
        }
        else
        {
            return variable_ <= other.variable_;
        }
    }

    bool operator >(const Var& other)
    {
        if (typeid(T) == typeid(std::string) || typeid(T) == typeid(char))
        {
            return variable_.length() > other.variable_.length();
        }
        else
        {
            return variable_ > other.variable_;
        }
    }
    bool operator >=(const Var& other)
    {
        if (typeid(T) == typeid(std::string) || typeid(T) == typeid(char))
        {
            return variable_.length() >= other.variable_.length();
        }
        else
        {
            return variable_ >= other.variable_;
        }
    }

    bool operator ==(const Var& other)
    {
        if (typeid(T) == typeid(std::string) || typeid(T) == typeid(char))
        {
            return variable_.length() == other.variable_.length();
        }
        else
        {
            return variable_ == other.variable_;
        }
    }

    bool operator !=(const Var& other)
    {
        if (typeid(T) == typeid(std::string) || typeid(T) == typeid(char))
        {
            return variable_.length() != other.variable_.length();
        }
        else
        {
            return variable_ != other.variable_;
        }
    }

private:
    T variable_;
};

int main()
{

}
